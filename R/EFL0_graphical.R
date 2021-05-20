#' @title EFL0_graphical
#' 
#' @description Elastic Fuse L0 (EFL0) Network.
#' 
#' @param data_current input dataset. Must include a column \code{weight} and a
#' column for each node in \code{vars}. 
#' @param vars vector of node names.
#' @param rlambda fraction of \code{lambda.max} to determine the smallest value 
#' for \code{lambda}. The default is 0.0001.
#' @param nlambda number of \code{lambda} values. Default is 50.
#' @param nlambda2 number of \code{lambda2} values. Default is 6.
#' @param beta_diff_max optional maximum tolerated change in beta coefficients 
#' from one time point to the next. This is used to determine the sequence of 
#' candidate \code{lamdbda2} values for the L2 penalty. Default is 1.
#' @param adaptive logical indicator of whether the overall dynamic model 
#' includes adaptive lasso weights for the L1 penalty. Default is \code{FALSE}.
#' @param adaptive_current logical indicator of whether to include adaptive 
#' lasso weights for the L1 penalty in this time point's model. Default is 
#' \code{FALSE}.
#' @param adaptive_current_wbeta vector of non-negative coefficient-specific 
#' adaptive lasso weights for the L1 penalty. Default is \code{NULL}.
#' @param adaptive_ridge logical indicator of whether the adaptive lasso weights
#' should be estimated with ridge regression. Default is \code{FALSE}, which 
#' means they will be estimated with OLS regression.
#' @param L0 logical indicator for whether to perform L0 penalty via 
#' hard-thresholding. Default is \code{TRUE}.
#' @param ncutB number of cuts to try for hard-thresholding. Default is 10.
#' @param EFL0 logical indicator of whether the overall dynamic model 
#' includes the two EFL0 penalties (L2 and L0). Default is \code{FALSE}.
#' @param EFL0_current logical indicator of whether to include the two EFL0 
#' penalties (L2 and L0) in this time point's model. Default is \code{FALSE}.
#' @param EFL0_current_ridgeC vector of coefficient-specific constant C values 
#' for the L2 penalty of EFL0. Default is\code{NULL}.
#' @param lambda2_seq_input vector of candidate \code{lambda2} values for the L2 
#' penalty of EFL0. Default is \code{NULL}.
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param foldid_subj an optional vector of values between 1 and \code{nfolds} 
#' specifying which fold each subject (or observation if there is only 1 
#' subject) is in.
#' @param precision_edges_true_current optional tibble of true edges for 
#' performance evaluation. Must include columnsm\code{row}, \code{col}, and
#' \code{edge}. Default is \code{NULL}.
#' @param seed randomization seed.
#' 
#' @return A list containing the following:
#' \item{coefs_adaptive_weights}{coefficients for adaptive lasso weights at the
#' next time point.}
#' \item{ridgeC_elastic_fuse}{final beta coefficients to be used in the EFL0
#' penalty at the next time point.}
#' \item{precision_edges}{a long dataset of precision matrix elements.}
#' \item{partial_edges}{a long dataset of partial correlation matrix elements.}
#' \item{edge_detection}{a matrix of tuning parameter values and performance 
#' metrics.}
#' \item{output_steps}{a long dataset of beta coefficients after lasso (L1), 
#' after soft-thresholding (L2), and after hard-thresholding (L0), when 
#' applicable.}
#' 
#' @export

################################################################################
# Program: EFL0_graphical.R
# Project: Smooth Time Varying Network
# Author: Erin McDonnell
# Date: 02-11-2020
#
# This function calls the EFL0() function to fit the p regression models 
# necessary to conduct glasso.  
# Depending on the options selected, we choose the optimal lambda parameters 
# from a grid, which is unique 
# for each of the p regressions.

EFL0_graphical = function(data_current, 
                          vars,
                          rlambda = 0.0001, 
                          nlambda = 50,
                          nlambda2 = 6,
                          beta_diff_max = 1,
                          adaptive = FALSE,
                          adaptive_current = FALSE,
                          adaptive_current_wbeta = NULL,
                          adaptive_ridge = FALSE,
                          L0 = FALSE, 
                          ncutB = 7,
                          EFL0 = FALSE,
                          EFL0_current = FALSE,
                          EFL0_current_ridgeC = NULL,
                          lambda2_seq_input = NULL,
                          nfolds = 5,
                          foldid_subj = NULL,
                          beta_min_threshold = 0,
                          precision_edges_true_current, 
                          seed = 123){
  set.seed(seed)
  
  p = length(vars)
  
  foldid = left_join(data_current, foldid_subj) %>%
    pull(foldid)
  
  weights_current = data_current$weight
  
  data_current = as.matrix(data_current %>% select(all_of(vars)))
  
  weights_current_scaled = weights_current * nrow(data_current) / sum(weights_current)
  
  weighted_mean = NULL
  weighted_sd = NULL
  data_centered = NULL
  data_standardized = NULL
  data_standardized_weighted = NULL
  
  for(j in 1:p){
    weighted_mean[j]  = sum(weights_current_scaled * data_current[,j]) / sum(weights_current_scaled)
    weighted_sd[j]    = sqrt(sum(weights_current_scaled * ((data_current[,j] - weighted_mean[j])^2)) / sum(weights_current_scaled))
    data_centered     = cbind(data_centered, data_current[,j] - weighted_mean[j])
    data_standardized = cbind(data_standardized, (data_current[,j] - weighted_mean[j]) / weighted_sd[j])
    data_standardized_weighted = cbind(data_standardized_weighted, data_standardized[,j] * sqrt(weights_current_scaled))
  }
  
  #########################################################################################
  # Set up beta weights for adaptive lasso
  
  if(adaptive_current == TRUE){
    wbeta = adaptive_current_wbeta
  } else{
    wbeta = matrix(1, nrow = p, ncol = p)
  }
  
  #########################################################################################
  # Set up ridgeC for EFL0 penalty term
  
  if(EFL0_current){
    ridgeC = EFL0_current_ridgeC
  } else{
    lambda2_seq = NULL
    ridgeC = matrix(0, nrow = p, ncol = p)
  }

  #########################################################################################
  # Loop through the p variables and run EFL0 with the same sequence of 50 lambdas for each variable, 
  # Sum over the cvm's across all nodes within a given lambda.
  # Choose the lambda with the minimum cvm.
  # NOTE: For models with L0 penalty, we have to run the EFL0 function for each lambda value separately. 
  # Otherwise we can't get the Beta0 values for each lambda...only the optimal one.
  
  # Shells for output
  
  lambda = NULL
  lambdamax = NULL
  lambda_is_lambdamin = NULL
  lambda_is_lambdamax = NULL
  lambda2 = NULL
  lambda2_is_lambda2min = NULL
  lambda2_is_lambda2max = NULL
  # EM 2021-05-12: Add tibbles for lambda sequences and chosen lambdas for each outcome node.
  lambda_seqs = NULL
  chosen_lambdas = NULL
  cv_betas = list()
  cv_int = NULL
  betas_adaptive_weights = list()
  model_specific_cvm = NULL
  # EM 2021-03-19: Add output_steps
  output_steps = NULL
  
  for(j in 1:p){
    
    wbeta_scaled = as.numeric(wbeta[j, -j] * (p-1) / sum(wbeta[j,-j]))

    # Regularized model
    
    # If penalty = "Lasso", then we just need to run EFL0 once.
    if(EFL0_current == FALSE){
      
      print("Running model with no EFL0 arguments")
      
      temp_model = EFL0(x         = data_standardized[, -j], 
                        y         = data_standardized[,  j],
                        weights   = weights_current,
                        nlambda   = nlambda,
                        rlambda   = rlambda,
                        wbeta     = wbeta_scaled,
                        iL0       = L0,
                        ncutB     = ncutB,
                        nfolds    = nfolds,
                        foldid    = foldid,
                        keep.beta = FALSE,
                        isd       = FALSE,
                        iysd      = FALSE,
                        print_test = FALSE)
      
      # If EFL0 = TRUE, then we need to identify lambda1 and lambda2 sequences.
    } else {
      
      # Select lambda2_seq based on lambda1_seq for this node
      lambda1max = maxLambdaLmC(X = data_standardized_weighted[, -j], 
                                y = data_standardized_weighted[,  j],
                                alpha = 1,
                                wbeta = wbeta_scaled,
                                N0 = nrow(data_standardized_weighted),
                                p = eval(p - 1))
      
      lambda1_seq = lambda1max * (rlambda)^(c(0:(nlambda-1))/(nlambda-1))
      #lambda2_seq = seq(lambda1max, rlambda, length.out = nlambda2)
      if(is.null(lambda2_seq_input)){
        lambda2_seq = seq(0, beta_diff_max * max(lambda1_seq) * max(wbeta_scaled) / 2, length.out = nlambda2)
      } else{
        lambda2_seq = lambda2_seq_input
      }
      
      print("Running model with EFL0 arguments")
      print(paste0("EFL0: ", EFL0))
      print(paste0("EFL0_current: ", EFL0_current))
      print("ridgeC_EFL0:")
      print(as.numeric(ridgeC[j, -j]))
      print("Lambda2 Sequence:")
      print(lambda2_seq)
      print("Adaptive lasso weights:")
      print(wbeta_scaled)
      print("ncutB:")
      print(ncutB)
      
      temp_model = EFL0(x            = data_standardized[, -j], 
                        y            = data_standardized[,  j],
                        weights      = weights_current,
                        #lambda       = lambda1_seq,
                        nlambda      = nlambda,
                        rlambda      = rlambda,
                        wbeta        = wbeta_scaled,
                        iL0          = L0,
                        ncutB        = ncutB,
                        iEFL0        = EFL0_current,
                        lambda2_EFL0 = lambda2_seq,
                        ridgeC_EFL0  = as.numeric(ridgeC[j, -j]),
                        nfolds       = nfolds,
                        foldid       = foldid,
                        keep.beta    = FALSE,
                        isd          = FALSE,
                        iysd         = FALSE,
                        print_test   = FALSE)
    }
    
    # EM 2021-03-19: Add a "steps" dataset to the output, which includes Lasso, Soft-thresholding, and 
    #                Hard-thresholding results from the EFL0 model    
    
    output_steps = bind_rows(output_steps, 
                             tibble(step = "Lasso",
                                    coef = temp_model$Beta, 
                                    outcome = vars[j],
                                    variable = c(vars[-j])))
    
    if(L0){
      
      if(EFL0_current){
      
        output_steps = bind_rows(output_steps,  
                                 tibble(step = "Soft",
                                        coef = temp_model$Beta_soft, 
                                        outcome = vars[j],
                                        variable = c(vars[-j])))
      }
        
      output_steps = bind_rows(output_steps,
                               tibble(step = "Hard",
                                      coef = temp_model$Beta0,
                                      outcome = vars[j],
                                      variable = c(vars[-j])))
    }
    # EM 2021-03-19: End update

    if(EFL0_current == TRUE){
      lambda[j] = temp_model$lambda1.opt
      lambda2[j] = temp_model$lambda2.opt
      lambda2_is_lambda2min[j] = as.numeric(lambda2[j] == min(lambda2_seq))
      lambda2_is_lambda2max[j] = as.numeric(lambda2[j] == max(lambda2_seq))
      cv_betas[[j]] = temp_model$Beta0
      model_specific_cvm = bind_rows(model_specific_cvm, tibble(outcome = paste0("node_", j), cvm = temp_model$fit0$cvm))
    } else if(L0 == TRUE){
      lambda[j] = temp_model$lambda.opt
      cv_betas[[j]] = temp_model$Beta0
      model_specific_cvm = bind_rows(model_specific_cvm, tibble(outcome = paste0("node_", j), cvm = temp_model$fit0$cvm))
    } else{
      lambda[j] = temp_model$lambda.min
      cv_betas[[j]] = temp_model$Beta
      model_specific_cvm = bind_rows(model_specific_cvm, tibble(outcome = paste0("node_", j), cvm = min(temp_model$fit$cvm)))
    }
    lambdamax[j] = max(temp_model$fit$lambda)
    lambda_is_lambdamin[j] = as.numeric(lambda[j] == min(temp_model$fit$lambda))
    lambda_is_lambdamax[j] = as.numeric(lambda[j] == max(temp_model$fit$lambda))
    cv_int[j] = temp_model$a
    
    # EM 2021-05-19: Use the "beta min condition" to further threshold estimates at s * sqrt(log(p/n))
    #s = sum(cv_betas[[j]] != 0)
    #cv_betas[[j]][which(abs(cv_betas[[j]]) < s * sqrt(log(p/nrow(data_standardized))))] = 0
    cv_betas[[j]][which(abs(cv_betas[[j]]) < beta_min_threshold)] = 0
    
    # EM 2021-05-12: Add a tibble of lambda sequences for each outcome node
    lambda_seqs = bind_rows(lambda_seqs, tibble(lambda_seq = temp_model$fit$lambda, tuning_parm = "lambda1", outcome_node = j))
    if(EFL0_current){
      lambda_seqs = bind_rows(lambda_seqs, tibble(lambda_seq = lambda2_seq, tuning_parm = "lambda2", outcome_node = j))
    }

    # Non-regularized model for constructing the next time point's adaptive weights
    if(adaptive == TRUE){
      if(adaptive_ridge == TRUE){
      temp_model_adaptive_weights = cv.glmnet(x = data_standardized[, -j], 
                                              y = data_standardized[,  j], 
                                              weights = weights_current, 
                                              alpha = 0,
                                              nlambda = nlambda,
                                              lambda.min.ratio = rlambda,
                                              nfolds = nfolds, 
                                              foldid = foldid)
                                    
      betas_adaptive_weights[[j]] = temp_model_adaptive_weights$glmnet.fit$beta[,which(temp_model_adaptive_weights$lambda == temp_model_adaptive_weights$lambda.min)]
      } else{
      temp_model_adaptive_weights = lm(data_standardized[,  j] ~ data_standardized[, -j], weights = weights_current)
      betas_adaptive_weights[[j]] = temp_model_adaptive_weights$coefficients[-1]
      }
    }
  }

  if(EFL0_current == TRUE){
    lambda2_range = max(lambda2) - min(lambda2)
    lambda2_num_zero = sum(lambda2 == 0)
  } else{
    lambda2 = NA
    lambda2_range = NA
    lambda2_num_zero = NA
  }
  lambda_range = max(lambda) - min(lambda)
  lambdamax_range = max(lambdamax) - min(lambdamax)
  lambda_num_zero = sum(lambda == 0)
  
  # EM 2021-05-12: Add a tibble of chosen lambdas for each outcome node.
  if(EFL0_current){
    chosen_lambdas = bind_rows(chosen_lambdas, 
                               tibble(outcome_node = 1:p,
                                      lambda1 = lambda, 
                                      lambda2 = lambda2,
                                      lambda_is_lambdamin = lambda_is_lambdamin,
                                      lambda_is_lambdamax = lambda_is_lambdamax,
                                      lambda2_is_lambda2min = lambda2_is_lambda2min,
                                      lambda2_is_lambda2max = lambda2_is_lambda2max))
  } else{
    chosen_lambdas = bind_rows(chosen_lambdas, 
                               tibble(outcome_node = 1:p,
                                      lambda1 = lambda, 
                                      lambda_is_lambdamin = lambda_is_lambdamin,
                                      lambda_is_lambdamax = lambda_is_lambdamax))
  }

  #########################################################################################
  # Now focus on the p EFL0 models we need (corresponding to our selected lambdas).
  # Store the betas from lambda = 0 for the adaptive lasso weights at the next time point
  # Get variance from each model
  # sigma^2 = 1/n * SSE
  
  coefs = NULL
  coefs_adaptive_weights = NULL
  ridgeC_elastic_fuse = NULL
  sigma2s = vector()

  for(j in 1:p){
    
    # Unstandardize
    reg_int = cv_int[j] * weighted_sd[j]
    reg_betas = cv_betas[[j]] * weighted_sd[j] / weighted_sd[-j]
    
    coefs = bind_rows(coefs, 
                      tibble(coef = reg_betas, 
                             variable = vars[-j]) %>% 
                        pivot_wider(names_from = variable, 
                                    values_from = coef))

    if(adaptive == TRUE){
      coefs_adaptive_weights = bind_rows(coefs_adaptive_weights, 
                                         tibble(coef = as.vector(betas_adaptive_weights[[j]]), 
                                                variable = vars[-j]) %>% 
                                           pivot_wider(names_from = variable, 
                                                       values_from = coef))
    }
    
    if(EFL0 == TRUE){
      ridgeC_elastic_fuse = bind_rows(ridgeC_elastic_fuse, 
                                      tibble(coef = as.vector(cv_betas[[j]]), 
                                             variable = vars[-j]) %>% 
                                        pivot_wider(names_from = variable, 
                                                    values_from = coef))
    }
    
    reg_sigma2 = sum(weights_current * (reg_int + data_centered[,-j] %*% reg_betas - data_centered[,j])^2) / sum(weights_current)
    sigma2s = c(sigma2s, reg_sigma2)
  }
  
  #########################################################################################
  # Set up a diagonal matrix of 1/sigma^2, and a matrix of regression coefficients
  # Use these matrices to calculate the Precision matrix
  
  if(adaptive == TRUE){
    coefs_adaptive_weights = coefs_adaptive_weights %>%
    select(all_of(vars))
  }
  
  if(EFL0 == TRUE){
    ridgeC_elastic_fuse = ridgeC_elastic_fuse %>%
      select(all_of(vars))
  }
  
  coefs = coefs %>%
    select(all_of(vars))

  coefs = as.matrix(coefs)
  
  diag(coefs) = 0
  if(all(coefs == 0)){
    print("Note: All coefficients are 0 across the p regressions")
  }
  
  D = diag(x = 1/sigma2s, nrow = length(sigma2s), ncol = length(sigma2s))
  
  precision = D %*% (diag(nrow = p, ncol = p) - coefs)

  #########################################################################################
  # Make the precision matrix symmetric
  # If at least one of the edges is 0, then they're both set to 0.
  # Otherwise they're both set to the mean.
  # Construct both precision and partial correlation matrices and edge tibbles
  
  precision_edges = NULL
  partial_edges = NULL
  
  for(j in 2:p){
    for(k in 1:(j-1)){
      j_k = precision[j, k]
      k_j = precision[k, j]
      precision[j, k] = precision[k, j] = ifelse(sign(j_k) == sign(k_j), sign(j_k) * sqrt(j_k * k_j), 0)
    }
  }
  
  partial = matrix(NA, p, p)
  
  for(j in 2:p){
    for(k in 1:(j-1)){
      partial[j, k] = partial [k, j] = -precision[j, k] / sqrt(precision[j,j] * precision[k,k])
    }
  }
  diag(partial) = 1
  
  # Generate column names
  colnames(precision) = seq(1, p)
  colnames(partial)   = seq(1, p)
  
  precision_edges = as_tibble(precision) %>% 
    mutate(row = seq(1, p)) %>%
    pivot_longer(names_to = "col",
                 values_to = "edge",
                 seq(1, p)) %>%
    mutate(col = as.numeric(col)) %>%
    select(row, col, edge)
  
  partial_edges = as_tibble(partial) %>% 
    mutate(row = seq(1, p)) %>%
    pivot_longer(names_to = "col",
                 values_to = "edge",
                 seq(1, p)) %>%
    mutate(col = as.numeric(col)) %>%
    select(row, col, edge)
  
  if(is.null(precision_edges_true_current)){
    return(list(coefs_adaptive_weights = coefs_adaptive_weights,
                ridgeC_elastic_fuse = ridgeC_elastic_fuse,
                precision_edges = precision_edges,
                partial_edges = partial_edges,
                # EM 2021-05-19: Add model-specific cvm
                model_specific_cvm = model_specific_cvm,
                # EM 2021-03-19: Add output_steps
                output_steps = output_steps,
                # EM 2021-05-12: Add node-level lambda information
                lambda_seqs = lambda_seqs,
                chosen_lambda = chosen_lambdas))
  } else{
    
    precision_edges_comparison = full_join(precision_edges_true_current %>% arrange(row, col) %>% rename(edge_true = edge),
                                           precision_edges              %>% arrange(row, col) %>% rename(edge_est  = edge),
                                           by = c("row", "col"))
    
    precision_edges_comparison_unique_edges = precision_edges_comparison %>%
      filter(row < col)

    edgetp = nrow(precision_edges_comparison_unique_edges %>% filter(edge_true != 0 & edge_est != 0)) 
    edgefn = nrow(precision_edges_comparison_unique_edges %>% filter(edge_true != 0 & edge_est == 0))
    edgefp = nrow(precision_edges_comparison_unique_edges %>% filter(edge_true == 0 & edge_est != 0))
    edgetn = nrow(precision_edges_comparison_unique_edges %>% filter(edge_true == 0 & edge_est == 0)) 
    sse    =  sum(precision_edges_comparison_unique_edges %>% mutate(squared_error = (edge_true - edge_est)^2) %>% pull(squared_error))

    return(list(coefs_adaptive_weights = coefs_adaptive_weights,
                ridgeC_elastic_fuse = ridgeC_elastic_fuse,
                precision_edges = precision_edges,
                partial_edges = partial_edges,
                edge_detection = tibble(tp     = edgetp, 
                                        fn     = edgefn, 
                                        fp     = edgefp, 
                                        tn     = edgetn,
                                        sse    = sse,
                                        cvm    = sum(model_specific_cvm$cvm)),
                # EM 2021-05-19: Add model-specific cvm
                model_specific_cvm = model_specific_cvm,
                # EM 2021-03-19: Add output_steps
                output_steps = output_steps,
                # EM 2021-05-12: Add node-level lambda information
                lambda_seqs = lambda_seqs,
                chosen_lambdas = chosen_lambdas))
  }
}