// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  int i, p=X.cols(), N=X.rows();
  Eigen::VectorXd mX(p), sdX(p);
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  return List::create(Named("x")=X, Named("sd")=sdX, Named("m")=mX);
}


/////////////////////////////////
/////   Linear Regression   /////
/////////////////////////////////


/*****  LM: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i;
  double LiMax=0.0, LiMaxi=0.0;
  
  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(y.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }
  
  LiMax=LiMax/N0/alpha;
  
  return(LiMax);
}


/*****  LM: Enet (L1+L2)  *****/
// EM 2020-03-04: Add "weights" option
// [[Rcpp::export]]
List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd weights, 
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta, Eigen::VectorXd wbeta2,
             int p, int N0, double thresh, int maxit, double thresh2){
  
  int  i, j, it=0, il, iadd, ia=0;
  double lambda2, lambda2i, zi, obj0, obj1, b0, db0, objQi, objQj, rss0, rss1, RSS0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd RSS=Eigen::VectorXd::Zero(nlambda), RSQ=Eigen::VectorXd::Zero(nlambda);
  double xr, dbMax;
  // EM 2020-03-04: Add v(N0), vX(p)
  Eigen::VectorXd mX(p), v(N0), vX(p), sdX(p), di(p);
  
  double a0=0.0, my=0.0;
  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);
  
  // EM 2020-03-04: Get WEIGHTED mean and WEIGHTED sd, standardize accordingly
  //                Multiply by sqrt(weight)
  //                Updates are similar to the standard() Fortran subfunction in the glmnet package
  
  weights = weights.array() * N0 / weights.sum();
  for(i=0;i<N0;++i){
    v(i) = sqrt(weights(i));
  }
  
  for (i=0;i<p;++i) {
    mX(i) = weights.dot(X.col(i)) / N0;
    X.col(i) = X.col(i).array()-mX(i);
    for(j=0;j<N0;++j){
      X.col(i)(j) = v(j) * X.col(i)(j);
    }
    vX(i) = X.col(i).dot(X.col(i)) / N0;
    sdX(i) = sqrt(vX(i));
    X.col(i) /= sdX(i);
  }
  my = weights.dot(y) / N0;
  y = (y.array() - my);
  for(i=0;i<N0;++i){
    y(i) = v(i) * y(i);
  }
  // EM 2020-03-04: End update
  
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }
  
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    
    it=0;
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          lambda2i=lambda2*wbeta2(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2i+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2i;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2i+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2i;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N0/2.0+objQj+objQi/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while
      
      iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array()/sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    a0=my;
    for(i=0;i<ia;i++){
      j=active(i);
      a0-=mX(j)*Beta(j,il);
    }
    a0S(il)=a0;
    
    if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
    // EM 2020-03-06: Add X, y, weights_scaled, and sdX to the return output
    return(List::create(Named("X")=X, Named("y")=y, Named("weights_scaled")=weights, Named("sdX")=sdX, 
                        Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, 
                        Named("a")=a0S, Named("flag")=flag, Named("rsq")=RSQ,
                        Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}


/*****  LM: Enet (L1+L2) cross-validation  *****/
// EM 2020-03-04: Add "weights" option
// [[Rcpp::export]]
List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd weights, 
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::VectorXd wbeta2,
               int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, Eigen::VectorXd weightsF, 
               int NF, double thresh2){
  
  int i, j, it=0, il, iadd, ia=0;
  double lambda2, lambda2i, zi, obj0, obj1, rss0, rss1, b0, db0, objQi, objQj, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda), BetaSTD=Eigen::MatrixXd::Zero(p,nlambda); // beta matrix for different lambdas
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXd RSSp(nlambda), RSS(nlambda), RSQ(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  double xr, dbMax;
  // EM 2020-03-04: Add v(N), vF(NF), ss(NF), vx(p)
  Eigen::VectorXd mX(p), v(N), vF(NF), ss(NF), vX(p), sdX(p), di(p);
  double a0=0.0, my=0.0;
  Eigen::MatrixXd predY=Eigen::MatrixXd::Zero(NF, nlambda);
  
  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);
  
  // EM 2020-03-04: Get WEIGHTED mean and WEIGHTED sd, standardize accordingly
  //                Multiply by sqrt(weight)
  //                Updates are similar to the standard() Fortran subfunction in the glmnet package
  
  weightsF = weightsF.array() * NF / weightsF.sum();
  for(i=0;i<NF;++i){
    vF(i) = sqrt(weightsF(i));
  }
  
  weights = weights.array() * N / weights.sum();
  for(i=0;i<N;++i){
    v(i) = sqrt(weights(i));
  }
  
  for (i=0;i<p;++i) {
    mX(i) = weights.dot(X.col(i)) / N;
    X.col(i) = X.col(i).array()-mX(i);
    for(j=0;j<N;++j){
      X.col(i)(j) = v(j) * X.col(i)(j);
    }
    vX(i) = X.col(i).dot(X.col(i)) / N;
    sdX(i) = sqrt(vX(i));
    X.col(i) = X.col(i) / sdX(i);
  }
  my = weights.dot(y) / N;
  y = (y.array() - my);
  for(i=0;i<N;++i){
    y(i) = v(i) * y(i);
  }
  // EM 2020-03-04: End update
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    
    it=0;
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          lambda2i=lambda2*wbeta2(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2i+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2i;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2i+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2)*lambda2i;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N/2.0+objQj+objQi/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while
      
      iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array()/sdX.array(); RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    a0=my; xbF.setZero(NF);
    for(i=0;i<ia;i++){
      j=active(i);
      xbF+=XF.col(j)*Beta(j,il);
      a0-=mX(j)*Beta(j,il);
    }
    xbF=xbF.array()+a0;
    
    predY.col(il)=xbF;
    
    // EM 2020-03-04: Make sure that the RSS is weighted, using the weights of the held-out fold.
    ss = Eigen::VectorXd::Zero(NF);
    
    for(i=0;i<NF;i++){
      ss(i) = vF(i) * (yF(i) - xbF(i));
    }
    RSSp(il)= ss.squaredNorm();
    // EM 2020-03-04: End update
    a0S(il)=a0;
    
    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
    return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, 
                        Named("flag")=flag, Named("predY")=predY, 
                        Named("a0S")=a0S, Named("RSS")=RSS, Named("rsq")=RSQ, 
                        Named("RSSp")=RSSp, Named("nlambda")=il));
}


// EM 2020-06-17: Add soft-thresholding function for Elastic Fuse L0 penalty
// EM 2021-04-14: Set dimension of vX to be p, not N
/*****  Used for soft-thresholding Beta (Elastic Fuse L0 penalty) *****/
// [[Rcpp::export]]
List softLmC(Eigen::VectorXd beta, Eigen::VectorXd betaSTD, 
             double lambda1, Eigen::VectorXd wbeta, double lambda2_EFL0, Eigen::VectorXd ridgeC_EFL0, 
             Eigen::MatrixXd X,  Eigen::VectorXd weights,  int N, int p) {
  int i, j;
  Eigen::VectorXd  S(p), betaSSTD(p), betaS(p), mX(p), vX(p), sdX(p), diff(p), left(p), right(p);
  Eigen::VectorXd v(N);
  double my;
  
  // Get WEIGHTED mean and WEIGHTED sd to standardize betaS accordingly.
  
  weights = weights.array() * N / weights.sum();
  for(i=0;i<N;++i){
    v(i) = sqrt(weights(i));
  }
  
  for (j=0;j<p;++j) {
    mX(j) = weights.dot(X.col(j)) / N;
    X.col(j) = X.col(j).array()-mX(j);
    for(i=0;i<N;++i){
      X.col(j)(i) = v(i) * X.col(j)(i);
    }
    vX(j) = X.col(j).dot(X.col(j)) / N;
    sdX(j) = sqrt(vX(j));
  }
  
  S.setZero(p); 
  
  for (j=0; j<p; ++j) {
    // Soft-thresholding operator
    
    // EM 2021-04-14: Check the pieces of the soft-thresholding operator
    diff(j)  = ridgeC_EFL0(j) - betaSTD(j);
    left(j)  = (lambda1 * wbeta(j)) / (2 * lambda2_EFL0);
    right(j) = std::abs(ridgeC_EFL0(j) - betaSTD(j));
    
    // EM 2021-03-11: slight modification to use absolute value instead of < - (should not change any results)
    if ((ridgeC_EFL0(j) - betaSTD(j)) > 0 & (lambda1 * wbeta(j)) / (2 * lambda2_EFL0) < std::abs(ridgeC_EFL0(j) - betaSTD(j))) {
      S(j) = ridgeC_EFL0(j) - betaSTD(j) - (lambda1 * wbeta(j)) / (2 * lambda2_EFL0);
    } else if ((ridgeC_EFL0(j) - betaSTD(j)) < 0 & (lambda1 * wbeta(j)) / (2 * lambda2_EFL0) < std::abs(ridgeC_EFL0(j) - betaSTD(j))){
      S(j) = ridgeC_EFL0(j) - betaSTD(j) + (lambda1 * wbeta(j)) / (2 * lambda2_EFL0);
    }
    // Final estimate is original beta + soft thresholding operator
    betaSSTD(j) = betaSTD(j) + S(j);
    
    // Get betas on original scale and calculate intercept and XB
    betaS(j) = betaSSTD(j) / sdX(j);
  }
  
  // EM 2021-03-11: Return S
  // EM 2021-04-14: Return all parameters
  return(List::create(Named("Beta")=betaS, Named("BetaSTD")=betaSSTD, Named("S")=S,
                      Named("i")=i, Named("j")=j, Named("mX")=mX, Named("sdX")=sdX,
                      Named("v")=v, Named("vX")=vX, Named("my")=my,
                      Named("diff")=diff, Named("left")=left, Named("right")=right));
}

/*****  Used for CV hard-thresholding cut Beta *****/
// EM 2020-03-12: Add "weights" option
// [[Rcpp::export]]
List cvHardLmC(Eigen::VectorXd beta, Eigen::VectorXd betaSTD, Eigen::VectorXd cut, Eigen::VectorXd wbeta,
               Eigen::MatrixXd X,  Eigen::VectorXd y,  Eigen::VectorXd weights,  int N, int p,
               Eigen::MatrixXd XF, Eigen::VectorXd yF, Eigen::VectorXd weightsF, int NF) {
  int i, j, nc=cut.size();
  Eigen::VectorXd betaC(p), mX(p);
  Eigen::VectorXd v(N), vF(NF), RSS, RSSp, xbF=Eigen::VectorXd::Zero(NF), xb=Eigen::VectorXd::Zero(N),
    ss=Eigen::VectorXd::Zero(N), ssp=Eigen::VectorXd::Zero(NF);
  double my, a0;
  
  // EM 2020-03-18: Scale weights for both the training and testing folds
  weightsF = weightsF.array() * NF / weightsF.sum();
  for(i=0;i<NF;++i){
    vF(i) = sqrt(weightsF(i));
  }
  
  weights = weights.array() * N / weights.sum();
  for(i=0;i<N;++i){
    v(i) = sqrt(weights(i));
  }
  // EM 2020-03-04: End update
  
  // EM 2020-03-18: Get WEIGHTED means of y and X
  for (i=0;i<p;++i) {
    mX(i) = weights.dot(X.col(i)) / N;
  }
  my = weights.dot(y) / N;
  // EM 2020-03-04: End update
  
  RSS.setZero(nc); RSSp.setZero(nc);
  
  for (i=0; i<nc; ++i) {
    
    betaC.setZero(p); a0=my;
    xb.setZero(N); xbF.setZero(NF);
    
    for (j=0; j<p; ++j) {
      if (std::abs(betaSTD(j)) > cut(i) || wbeta(j)==0.0) {
        betaC(j)=beta(j);
        a0-=mX(j)*betaC(j);
        
        xb+=X.col(j)*betaC(j);
        xbF+=XF.col(j)*betaC(j);
      }
    }
    xb=xb.array()+a0;
    xbF=xbF.array()+a0;
    
    // EM 2020-03-18: Weight each residual by its square-root weight before calculating RSS and RSSp
    ss.setZero(N);
    for(j=0;j<N;j++){
      ss(j) = v(j) * (y(j) - xb(j));
    }
    RSS(i) = ss.squaredNorm();
    
    ssp.setZero(NF);
    for(j=0;j<NF;j++){
      ssp(j) = vF(j) * (yF(j) - xbF(j));
    }
    RSSp(i) = ssp.squaredNorm();
    // EM 2020-03-18: End update
    
  }
  
  return(List::create(Named("RSS")=RSS, Named("RSSp")=RSSp));
}

/*****  Used for CV trimming of number  *****/
// [[Rcpp::export]]
List cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
               Eigen::MatrixXd X,  Eigen::VectorXd y,  Eigen::VectorXd weights,  int N, int p,
               Eigen::MatrixXd XF, Eigen::VectorXd yF, Eigen::VectorXd weightsF, int NF) {
  int i, j, k;
  Eigen::VectorXd mX(p);
  Eigen::VectorXd v(N), vF(NF), RSS, RSSp, xbF=Eigen::VectorXd::Zero(NF), xb=Eigen::VectorXd::Zero(N),
    ss=Eigen::VectorXd::Zero(N), ssp=Eigen::VectorXd::Zero(NF);
  double my, a0;
  
  // EM 2020-03-18: Scale weights for both the training and testing folds
  weightsF = weightsF.array() * NF / weightsF.sum();
  for(i=0;i<NF;++i){
    vF(i) = sqrt(weightsF(i));
  }
  
  weights = weights.array() * N / weights.sum();
  for(i=0;i<N;++i){
    v(i) = sqrt(weights(i));
  }
  // EM 2020-03-18: End update  
  
  // EM 2020-03-18: Get WEIGHTED means of y and X
  for (i=0;i<p;++i) {
    mX(i) = weights.dot(X.col(i)) / N;
  }
  my = weights.dot(y) / N;
  // EM 2020-03-18: End update
  
  RSS.setZero(nn2); RSSp.setZero(nn2);
  
  if (nn==0) {
    
    a0=my;
    xb.setZero(N); xbF.setZero(NF);
    
    xb=xb.array()+a0;
    xbF=xbF.array()+a0;
    
    // EM 2020-03-18: Weight each residual by its square-root weight before calculating RSS and RSSp
    ss.setZero(N);
    for(k=0;k<N;k++){
      ss(k) = v(k) * (y(k) - xb(k));
    }
    RSS(0) = ss.squaredNorm();
    
    ssp.setZero(NF);
    for(k=0;k<NF;k++){
      ssp(k) = vF(k) * (yF(k) - xbF(k));
    }
    RSSp(0) = ssp.squaredNorm();
    // EM 2020-03-18: End update
    
  } else {
    
    a0=my;
    xb.setZero(N); xbF.setZero(NF);
    
    for (i=0; i<nn; ++i) {
      j=loco(i);
      a0-=mX(j)*beta(i);
      
      xb+=X.col(j)*beta(i);
      xbF+=XF.col(j)*beta(i);
      
      xb=xb.array()+a0;
      xbF=xbF.array()+a0;
      
      // EM 2020-03-18: Weight each residual by its square-root weight before calculating RSS and RSSp
      ss.setZero(N);
      for(k=0;k<N;k++){
        ss(k) = v(k) * (y(k) - xb(k));
      }
      RSS(i) = ss.squaredNorm();
      
      ssp.setZero(NF);
      for(k=0;k<NF;k++){
        ssp(k) = vF(k) * (yF(k) - xbF(k));
      }
      RSSp(i) = ssp.squaredNorm();
      // EM 2020-03-18: End update
      
      xb=xb.array()-a0;
      xbF=xbF.array()-a0;
    }
    
  }
  
  if (nn2>nn && nn==0) {
    for(i=nn;i<nn2;i++){RSS(i)=RSS(0);RSSp(i)=RSSp(0);}
  }
  if (nn2>nn && nn>0) {
    for(i=nn;i<nn2;i++){RSS(i)=RSS(nn-1);RSSp(i)=RSSp(nn-1);}
  }
  
  return(List::create(Named("RSS")=RSS, Named("RSSp")=RSSp));
}


