#' @title EFL0
#' 
#' @description Elastic Fuse L0 (EFL0) regression model.
#' 
#' @param x input matrix. Each row is an observation vector.
#' @param y response vector.
#' @param weights vector of observation-level weights.
#' @param lambda a use-supplied vector of candidate lambda values. If 
#' \code{NULL} (default), a sequence of \code{lambda} is generated based on 
#' \code{nlambda} and \code{rlambda}. Supplying a value of \code{lambda} 
#' overrides this.
#' @param nlambda number of \code{lambda} values. Default is 50.
#' @param rlambda fraction of \code{lambda.max} to determine the smallest value 
#' for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of 
#' observations is larger than or equal to the number of variables; otherwise, 
#' \code{rlambda = 0.01}.
#' @param wbeta vector of non-negative coefficient-specific adaptive lasso 
#' weights for the L1 penalty. Default is 1 for all coefficients.
#' @param nfolds number of folds. With \code{nfolds = 1} and 
#' \code{foldid = NULL} by default, cross-validation is not performed. 
#' Specifying \code{foldid} overrides \code{nfolds}.
#' @param foldid an optional vector of values between 1 and \code{nfolds} 
#' specifying which fold each observation is in.
#' @param iEFL0 logical indicator of whether to include the two EFL0 penalties 
#' (L2 and L0). Default is \code{FALSE}.
#' @param lambda2_EFL0 vector of candidate lambda2 values for the L2 penalty of 
#' EFL0. Default is \code{NULL}.
#' @param ridgeC_EFL0 vector of coefficient-specific constant C values for the 
#' L2 penalty of EFL0. Default is\code{NULL}.
#' @param iL0 logical flag for performing L0-norm via hard-thresholding. Default 
#' is \code{TRUE}. 
#' @param ncutB number of cuts to try for hard-thresholding. Default is 10.
#' @param isd logical flag for outputting standardized coefficients. \code{x} is
#' always standardized prior to fitting the model. Default is 
#' \code{isd = FALSE}, returning \eqn{\beta} on the original scale.
#' @param iysd logical flag for standardizing \code{y} prior to computation. The
#' returning coefficients are always based the original \code{y} 
#' (unstandardized). Default is \code{isd = FALSE}.
#' @param keep.beta logical flag for returning estimates for all \code{lambda} 
#' values. Default is \code{FALSE}.
#' @param thresh convergence threshold for coordinate descent. Default value is 
#' \code{1E-6}.
#' @param maxit Maximum number of iterations for coordinate descent. Default is 
#' \code{10^5}.
#' @param print_test logical indicator for whether to print test output. Default
#' is \code{FALSE}.
#' 
#' @return A list containing the following:
#' \item{a}{the intercept}
#' \item{Beta}{a sparse Matrix of coefficients after L1 penalization.}
#' \item{Beta_soft}{coefficients after additionally performing soft-thresholding
#' for \code{iEFL0 = TRUE}}
#' \item{Beta0}{coefficients after additionally performing L0-norm for 
#' \code{iL0 = TRUE}}
#' \item{fit}{a data.frame containing \code{lambda} and the number of non-zero 
#' coefficients \code{nzero}. 
#' With cross-validation, additional results are reported, such as average 
#' cross-validation partial likelihood \code{cvm} and its standard error 
#' \code{cvse}, and \code{index} with `*' indicating the minimum \code{cvm}.}
#' \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and 
#' \code{nzero} based on \code{iL0 = TRUE}. \code{cvm} in \code{fit0} may be 
#' different from \code{cvm} in \code{fit}, because of hard-thresholding.}
#' \item{lambda.opt}{value of \code{lambda} based on \code{iL0 = TRUE}.}
#' \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means 
#' converged.}
#' 
#' @export

#####################################################################################
# Code modified from the LmL0 function of the APML0 package (Li et al. 2018) 
# Program: EFL0.R
# Project: EFL0Network
# Author: Erin McDonnell
# Date created: 3/4/2020
# Date modified: 3/2/2021
#####################################################################################

# EM 2020-03-04: Add "weights" argument for observation weights
# EM 2020-05-26: Add "ridgeC" argument for ridge penalty fusing betas to fixed constants
# EM 2020-06-17: Add "iEFL0", "lambda2_EFL0", and "ridgeC_EFL0" argument for elastic fuse L0 penalty
# EM 2021-03-02: Remove the following arguments:
#                - family (will always be "gaussian")
#                - penalty (will always be "Lasso")
#                - Omega
#                - sgn
#                - lambda2
#                - ridgeC
#                - ill (will always be FALSE)
#                - ifast (will always be FALSE)
#                - icutB (will always be TRUE: perform **hard-thresholding** to the coefficients 

EFL0 = function(x, y, weights, lambda=NULL, nlambda=50, 
                rlambda=NULL, wbeta=rep(1,ncol(x)),
                nfolds=1, foldid=NULL,  
                iEFL0=FALSE, lambda2_EFL0=NULL, ridgeC_EFL0=NULL,
                iL0=TRUE, ncutB=10, isd=FALSE, iysd=FALSE, 
                keep.beta=FALSE, thresh=1e-6, maxit=1e+5,
                print_test = TRUE) {

  if(iEFL0){
    if(!iL0){
      stop("iEFL0 can only be TRUE if iL0 is also true")
    }
    if(is.null(lambda2_EFL0)){
      stop("iEFL0 set to TRUE but no value provided for lambda2_EFL0")
    }
    if(is.null(ridgeC_EFL0)){
      stop("iEFL0 set to TRUE but no value provided for ridgeC_EFL0")
    }
  } else{
    if(!is.null(lambda2_EFL0)){
      stop("lambda2_EFL0 provided but iEFL0 set to FALSE.")
    }
    if(!is.null(ridgeC_EFL0)){
      stop("ridgeC_EFL0 provided but iEFL0 set to FALSE.")
    }
  }
  
  # EM 2020-03-04: Add weight of 1 for each observation if no weights are supplied
  if (is.null(weights)){
    weights = rep(1, nrow(x))
  } else{
    if (length(weights) != nrow(x)) 
      stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of x (", nrow(x), ")", sep = ""))
  }
  
  wbeta=abs(wbeta)
  
  N0=nrow(x); p=ncol(x)
  
  sdy=1.0
  if (iysd) {
    sdy= sqrt(sum(weights * (y - weighted.mean(y, w=weights))^2) / sum(weights))
    y=y/sdy
  }
  
  ### Adaptive
  aPen=ifelse(all(wbeta>0), TRUE, FALSE)
  wbeta2=ifelse(wbeta==0,0,1)
  
  ### Lambda path
  if (all(wbeta==0)) {
    lambda=0.0
  }
  
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }
  
  #####  Run  #####
  # EM 2020-03-04: Add "weights" argument for penalty = "Net"
  # EM 2020-05-26: Add a "FuseEnetLmC" function call
  # EM 2020-05-30: Modify "FuseEnetLmC" function call
  out=EnetLmC(x, y, weights, alpha=1, lambda, nlambda, ilambda, wbeta, wbeta2, p, N0, thresh, maxit, 1e-5)

  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]
  
  out$Beta[is.na(out$Beta)]=out$BetaSTD[is.na(out$Beta)]
  out$Beta=Matrix(out$Beta[, 1:nlambdai,drop=FALSE]*sdy, sparse=TRUE)
  out$BetaSTD=Matrix(out$BetaSTD[, 1:nlambdai,drop=FALSE]*sdy, sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]
  out$a=out$a[1:nlambdai]
  
  
  if (nfolds==1 & is.null(foldid)) {
    
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    if (!isd) {
      return(list(X=out$X, y=out$y, weights_scaled=out$weights_scaled, sdX=out$sdX, sdy=sdy, a=out$a, Beta=out$Beta, fit=fit, flag=out$flag))
    } else {
      return(list(X=out$X, y=out$y, weights_scaled=out$weights_scaled, sdX=out$sdX, sdy=sdy, a=out$a, Beta=out$BetaSTD, fit=fit, flag=out$flag))
    }
    
  } else {
    
    ########################################
    #####  Cross-validation estimates  #####
    
    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    # EM 2020-03-18: We need both the number of observations in each fold AND the sums of the weights in each fold.
    Nfveci = as.vector(tapply(rep(1, N0), foldid, sum))
    sumweightsi = as.vector(tapply(weights, foldid, sum))
    
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      
      outi[[i]]=cvEnetLmC(x[!temid, ,drop=FALSE], y[!temid], weights[!temid], alpha=1, lambdai, nlambdai, wbeta, wbeta2, N0i[i],p, thresh, maxit, x[temid, ,drop=FALSE], y[temid], weights[temid], Nf[i], 1e-5)
      
      outi[[i]]$Beta[is.na(outi[[i]]$Beta)]=outi[[i]]$BetaSTD[is.na(outi[[i]]$Beta)]
      
      
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda] ## for ith fold

    }
    
    # EM 2020-03-18: Swap weighti with Nfveci and sumweightsi as appropriate        
    cvRSS=cvRSS[, 1:nlambdai,drop=FALSE]
    cvraw=cvRSS/Nfveci; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=sumweightsi, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=sumweightsi, na.rm=TRUE)/(nfoldi-1))
    
    
    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)
    
    
    if (!iL0) {
      if (!keep.beta) {
        if (!isd) {
          return(list(a=out$a[indexi], Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD[, indexi], fit=temCV, lambda.min=lambdai[indexi], flag=out$flag))
        }
        
      } else {
        if (!isd) {
          return(list(a=out$a, Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD, fit=temCV, lambda.min=lambdai[indexi], flag=out$flag))
        }
      }
    }
    
    
    
    # EM 2020-06-17: Add soft-thresholding step for elastic fuse L0 penalty
    ##########################################################
    #####  Elastic Fuse penalty: Soft-thresholding step  #####
    
    nlambda2 = ifelse(iEFL0, length(lambda2_EFL0), 1)
    # EM 2020-07-15: Make cv.min a matrix instead of a vector (grid search for both lambda1 and lambda2)
    cv.min = matrix(NA, nrow = nlambdai, ncol = nlambda2)
    cvm = list()
    
    # EM 2021-03-03: Generate soft-thresholded betas for every lambda1/lambda2 combination on the full dataset
    #                Create an out_soft object: A list of length nlambda2, each item is a list of 2 p-1 x nlambda1 matrices (1 for Beta, 1 for BetaSTD).
    out_soft = list()
    if(iEFL0){
      
      for(j in 1:nlambda2){
        
        Betaj_soft    = matrix(NA, nrow = p, ncol = nlambda)
        BetaSTDj_soft = matrix(NA, nrow = p, ncol = nlambda)
        
        for(k in 1:nlambdai){
          
          # EM 2021-03-11: Check that we have no negative thresholds
          if(any((lambdai[k] * wbeta / (2 * lambda2_EFL0[j])) < 0)){
            stop(paste("soft-threshold tuning parameters < 0 at (lambda1, lambda2) = ", lambdai[k], lambda2_EFL0[j]))
          }
          
          Betajk    = out$Beta[, k]
          BetaSTDjk = out$BetaSTD[, k]
          
          outjk_soft = softLmC(Betajk, BetaSTDjk, lambdai[k], wbeta, lambda2_EFL0[j], 
                               ridgeC_EFL0, x, weights, N0, p)
          
          Betaj_soft[,k]    = outjk_soft$Beta
          BetaSTDj_soft[,k] = outjk_soft$BetaSTD
        }
        
        # EM 2021-03-17: if lambda2_EFL0[j] = 0 and Betaj_soft != Betajk, then stop.
        if(lambda2_EFL0[j] == 0 & any(abs(Betaj_soft - out$Beta) > 0.0000001)){
          print("Betaj_soft")
          print(Betaj_soft)
          print("out$Beta")
          print(out$Beta)
          stop("Lambda2 = 0 but betas differ after soft-thresholding")
        }
        # EM 2021-03-17: End update.
        
        out_soft[[j]] = list(Beta = Betaj_soft, BetaSTD = BetaSTDj_soft)
      }
    } else{
      out_soft[[1]] = list(Beta = out$Beta, BetaSTD = out$BetaSTD)
    }
    # EM 2021-03-03: End Update.
    
    
    for(j in 1:nlambda2){
      
      if(print_test == TRUE){ 
        print(paste("j=", j)) 
      }
      
      #print(paste("j: ", j))
      
      if(iEFL0){
        
        outi_soft = list()
        
        for (i in 1:nfolds) {
          
          if(print_test == TRUE & j == 1){
            print(paste("i=", i))
          }
          
          temid = (foldid==i)
          
          if(print_test == TRUE & j == 1){
            print("outi[[i]]$Beta:")
            print(outi[[i]]$Beta)
            print("outi[[i]]$BetaSTD:")
            print(outi[[i]]$BetaSTD)
          }
          
          Betai_soft    = matrix(NA, nrow = p, ncol = nlambda)
          BetaSTDi_soft = matrix(NA, nrow = p, ncol = nlambda)
          
          for(k in 1:nlambdai){
            
            Betaij    = outi[[i]]$Beta[, k, drop=FALSE]
            BetaSTDij = outi[[i]]$BetaSTD[, k, drop=FALSE]
            
            if(print_test == TRUE & j == 1 & k >= 38){
              print(paste("j=1, i=", i, ", k=", k, ". Entering soft-thresholding."))
              print("Betaij")
              print(Betaij)
              print("BetaSTDij")
              print(BetaSTDij)
              print("lambdai[k]")
              print(lambdai[k])
              print("wbeta")
              print(wbeta)
              print("lambda2_EFL0[j]")
              print(lambda2_EFL0[j])
              print("ridgeC_EFL0")
              print(ridgeC_EFL0)
              print("x[!temid,,drop=FALSE]")
              print(head(x[!temid,,drop=FALSE]))
              print("weights[!temid]")
              print(weights[!temid])
              print("N0i[i]")
              print(N0i[i])
              print("p")
              print(p)
            }
              
            outij_soft = softLmC(Betaij, BetaSTDij, lambdai[k], wbeta, lambda2_EFL0[j], ridgeC_EFL0,
                                 x[!temid,,drop=FALSE], weights[!temid], N0i[i], p)
            
            if(print_test == TRUE & j == 1 & k >= 38){
              print("Soft-thresholding results:")
              print("outij_soft$Beta")
              print(outij_soft$Beta)
              print("outij_soft$BetaSTD")
              print(outij_soft$BetaSTD)
              print("outij_soft$S")
              print(outij_soft$S)
              print("i")
              print(outij_soft$i)
              print("j")
              print(outij_soft$j)
              print("mX")
              print(outij_soft$mX)
              print("sdX")
              print(outij_soft$sdX)
              print("v")
              print(outij_soft$v)
              print("vX")
              print(outij_soft$vX)
              print("my")
              print(outij_soft$my)
              print("diff")
              print(outij_soft$diff)
              print("left")
              print(outij_soft$left)
              print("right")
              print(outij_soft$right)
            }
          
            Betai_soft[,k] = outij_soft$Beta
            BetaSTDi_soft[,k] = outij_soft$BetaSTD
          }
            
          # EM 2021-03-17: if lambda2_EFL0[j] = 0 and Betaj_soft != Betajk, then stop.
          if(lambda2_EFL0[j] == 0 & any(abs(Betai_soft - outi[[i]]$Beta) > 0.0000001)){
            print("Betai_soft")
            print(Betai_soft)
            print("outi[[i]]$Beta")
            print(outi[[i]]$Beta)
            stop("Lambda2 = 0 but betas differ after soft-thresholding during cross-validation")
          }
          # EM 2021-03-17: End update.
            
          outi_soft[[i]] = list(Beta = Betai_soft, BetaSTD = BetaSTDi_soft)
        }
        # EM 2021-03-03: Use outi_soft in place of outi_new (just renamed)
      } else{
        outi_soft = outi
      }
      
      if(print_test == TRUE & j == 1){ 
        print("outi_soft") 
        print(outi_soft)
      }

      #####################################
      #####  Cross-validation for L0  #####

      # icutB = TRUE
      
      ###  cutoff  ###
      il0=1
      
      # EM 2020-07-15: Defined cv.min as a matrix outside the j loop.
      cvm[[j]]=list();
      repeat {
        
        if(print_test == TRUE & j == 1){ 
          print(paste("il0=", il0)) 
        }
        
        # EM 2021-03-03: Use the soft-thresholded betas from full dataset to determine cut points
        BetaSTD=out_soft[[j]]$BetaSTD[,il0]
        cut0=sort(unique(c(0,abs(BetaSTD))),decreasing=FALSE); cuti=NULL
        for (i in 1:ncutB) {
          cuti=c(cuti, diff(cut0)/ncutB*i+cut0[-length(cut0)])
        }
        cut0=sort(unique(c(cut0,cuti)))
        
        # EM 2021-03-10: Print cut thresholds for quality check
        #print("BetaSTD:")
        #print(BetaSTD)
        #print("Cut thresholds:")
        #print(cut0)
        
        # EM 2021-03-03: Use outi_soft in place of outi_new (just renamed)
        Betai=matrix(sapply(outi_soft, function(x){x$Beta[, il0,drop=FALSE]}), nrow=p)
        BetaSTDi=matrix(sapply(outi_soft, function(x){x$BetaSTD[, il0,drop=FALSE]}), nrow=p)
        
        cvRSS=matrix(NA, nrow=nfolds, ncol=length(cut0)); i=1
        for (i in 1:nfolds) {
          
          if(print_test == TRUE & j == 1 & il0 %in% c(1,nlambdai)){ 
            print(paste("i=", i)) 
          }
          
          temid=foldid==i
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          
          if(print_test == TRUE & j == 1 & il0 %in% c(1,nlambdai)){
            print("Betaj:")
            print(Betaj)
            print("BetaSTDj:")
            print(BetaSTDj)
          }
          
          # EM 2020-03-04: Add 'weights' option to cvHardLmC function call
          outij=cvHardLmC(Betaj, BetaSTDj, cut0, wbeta, 
                          x[!temid,,drop=FALSE], y[!temid], weights[!temid], N0i[i], p, 
                          x[ temid,,drop=FALSE], y[ temid], weights[ temid], Nf[i])
          
          cvRSS[i,]=outij$RSSp
        }
        
        if(print_test == TRUE & j == 1 & il0 %in% c(1,nlambdai)){
          print("cvRSS:")
          print(cvRSS)
        }
        
        # EM 2020-03-18: Swap weighti with Nfveci and sumweightsi as appropriate        
        cvraw=cvRSS/Nfveci; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
        cvm[[j]][[il0]]=apply(cvraw, 2, weighted.mean, w=sumweightsi, na.rm=TRUE)
        cv.min[il0,j]=min(cvm[[j]][[il0]])
        
        # EM 2021-03-10: Remove il1=which.min(cv.min[,j]) (il1 is not used anywhere else)
        if (nlambdai==1) break
        il0=il0+1
        if (il0>nlambdai) break
      }
      # End icutB = TRUE
    } # EM 2020-07-16: End loop through j
      
    # EM 2020-07-16: Pick best (lambda1, lambda2) combo and soft/hard threshold final betas accordingly.
    # EM 2020-07-22: If there are multiple (lambda1, lambda2) combos with the same minimum, pick the first combo.
    
    grid.min = as.matrix(which(cv.min == min(cv.min, na.rm = TRUE), arr.ind = TRUE))
    il01 = grid.min[1,1] # index of best lambda1 (first row means the first eligible lambda1,lambda2 combo)
    il02 = grid.min[1,2] # index of best lambda2 (first row means the first eligible lambda1,lambda2 combo)
    if(!iEFL0){
      print("lambda1,lambda2 indices:")
      print(il01)
      print(il02)
    }

    #print("cv.min matrix:")
    #print(cv.min)
    #print("which(cv.min == min(cv.min)):")
    #print(grid.min)
    #print(paste("il01: ", il01))
    #print(paste("il02: ", il02))
    nzero0=which.min(cvm[[il02]][[il01]])
    
    # 2021-03-03 EM: Save the soft estimates from the chosen lambda1/lambda2 combination.
    #                Perform final hard-thresholding on these betas.
    
    Beta_soft=out_soft[[il02]]$Beta[,il01]
    BetaSTD_soft=out_soft[[il02]]$BetaSTD[,il01]

    Beta0 = Beta_soft
    BetaSTD0 = BetaSTD_soft
    
    # icutB = TRUE
    
    cut0=sort(unique(c(0,abs(BetaSTD_soft))),decreasing=FALSE); cuti=NULL
    for (i in 1:ncutB) {
      cuti=c(cuti, diff(cut0)/ncutB*i+cut0[-length(cut0)])
    }
    cut0=sort(unique(c(cut0,cuti)))
    
    Beta0[abs(BetaSTD_soft)<=cut0[nzero0] & wbeta>0.0]=0.0
    BetaSTD0[abs(BetaSTD_soft)<=cut0[nzero0] & wbeta>0.0]=0.0
    
    # End icutB = TRUE
    
    # EM 2020-07-16: Final fit0 statistics
    # EM 2021-02-25: Add cvm_grid (cross-validation errors for all lambda1, lambda2 combos) to the output
    # EM 2021-03-11: Only ouput the lasso results that correspond to the OPTIMAL lamba1 (il01), not the minimum lambda1 from before soft/hard-thresholding (indexi).
    if(iEFL0){
      temCV0=data.frame(lambda1=lambdai[il01], lambda2=lambda2_EFL0[il02], cvm=cv.min[il01,il02], nzero=sum(Beta0!=0))
      
      if (!keep.beta) {
        
        if (!isd) {
          return(list(a=out$a[il01], Beta=out$Beta[, il01], Beta_soft=Beta_soft, Beta0=Beta0, fit=temCV, fit0=temCV0, 
                      lambda1.opt=lambdai[il01], lambda2.opt = lambda2_EFL0[il02], flag=out$flag, cvm_grid=cv.min))
        } else {
          return(list(Beta=out$BetaSTD[, il01], Beta_soft=BetaSTD_soft, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, 
                      lambda1.opt=lambdai[il01], lambda2.opt = lambda2_EFL0[il02], flag=out$flag, cvm_grid=cv.min))
        }
        
      } else {
        if (!isd) {
          return(list(a=out$a, Beta=out$Beta, Beta_soft=Beta_soft, Beta0=Beta0, fit=temCV, fit0=temCV0, 
                      lambda1.opt=lambdai[il0], lambda2.opt = lambda2_EFL0[il02], flag=out$flag, cvm_grid=cv.min))
        } else {
          return(list(Beta=out$BetaSTD, Beta_soft=BetaSTD_soft, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, 
                      lambda1.opt=lambdai[il0], lambda2.opt = lambda2_EFL0[il02], flag=out$flag, cvm_grid=cv.min))
        }
      }
    } else{
      temCV0=data.frame(lambda=lambdai[il01], cvm=cv.min[il01,il02], nzero=sum(Beta0!=0))
      
      if (!keep.beta) {
        
        if (!isd) {
          return(list(a=out$a[il01], Beta=out$Beta[, il01], Beta_soft=Beta_soft, Beta0=Beta0, fit=temCV, fit0=temCV0, 
                      lambda.opt=lambdai[il01], flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD[, il01], Beta_soft=BetaSTD_soft, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, 
                      lambda.opt=lambdai[il01], flag=out$flag))
        }
        
      } else {
        if (!isd) {
          return(list(a=out$a, Beta=out$Beta, Beta_soft=Beta_soft, Beta0=Beta0, fit=temCV, fit0=temCV0, 
                      lambda.opt=lambdai[il0], flag=out$flag))
        } else {
          return(list(Beta=out$BetaSTD, Beta_soft=BetaSTD_soft, Beta0=BetaSTD0, fit=temCV, fit0=temCV0, 
                      lambda.opt=lambdai[il0], flag=out$flag))
        }
      }
    }
  }
}




