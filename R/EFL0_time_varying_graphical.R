#' @title EFL0_time_varying_graphical
#' 
#' @description Dynamic Elastic Fuse L0 (EFL0) Network.
#' 
#' @param data input dataset. Must include columns \code{id}, \code{s_il}, and
#'  a column for each node in \code{vars}. 
#' @param vars vector of node names.
#' @param visits_of_interest vector of times at which to estimate the network.
#' @param visits_adaptive order in which to estimate the time-specific networks.
#' @param precision_edges_true optional tibble of true edges for performance 
#' evaluation. Must include columns \code{s_il}, \code{row}, \code{col}, 
#' \code{edge}. Default is \code{NULL}.
#' @param kernel logical indicator of whether to include observation-level 
#' kernel weights. Default is \code{FALSE}.
#' @param nh number of candidate bandwidths h to try for kernel weights. 
#' Default is 9.
#' @param h_power_low minimum power for constructing bandwidth sequence. 
#' Default is -0.5.
#' @param h_power_high maximum power for constructing bandwidth sequence. 
#' Default is 1.5.
#' @param adaptive logical indicator of whether to include adaptive lasso 
#' weights for the L1 penalty. Default is \code{FALSE}.
#' @param adaptive_ridge logical indicator of whether the adaptive lasso weights
#' should be estimated with ridge regression. Default is \code{FALSE}, which 
#' means they will be estimated with OLS regression.
#' @param L0 logical indicator of whether to include an L0 penalty. Default is 
#' \code{FALSE}.
#' @param EFL0 logical indicator of whether to include the Elastic Fuse 
#' penalties (L2 and L0). Default is \code{FALSE}.
#' @param beta_diff_max optional maximum tolerated change in beta coefficients 
#' from one time point to the next. This is used to determine the sequence of 
#' candidate \code{lamdbda2} values for the L2 penalty. Default is 1.
#' @param lambda2_seq_input user-inputted sequence of candidate \code{lambda2} 
#' values. Default is \code{NULL}.
#' @param nfolds number of folds. Default is 5.
#' @param foldid_subj an optional vector of values between 1 and \code{nfolds} 
#' specifying which fold each subject (or observation if there is only 1 
#' subject) is in.
#' @param seed randomization seed.
#' 
#' @return A list containing the following:
#' \item{visit_weights}{a dataset of observation-level kernel weights.}
#' \item{precision_edges}{a long dataset of precision matrix elements.}
#' \item{partial_edges}{a long dataset of partial correlation matrix elements.}
#' \item{edge_detection}{a matrix of tuning parameter values and performance 
#' metrics.}
#' \item{foldid}{the foldids that were used to perform cross-validation.}
#' \item{output_steps}{a long dataset of beta coefficients after lasso (L1), 
#' after soft-thresholding (L2), and after hard-thresholding (L0), when 
#' applicable.}
#' \item{h_selection}{a dataset of candidate bandwidths and resulting 
#' cross-validation errors.}
#' 
#' @export


#########################################################################################
# Program: EFL0_time_varying_graphical.R
# Project: Smooth Time Varying Network
# Author: Erin McDonnell
# Date: 02-11-2020
#
# This function calls the EFL0_graphical() function for each time point of interest,
# storing the previous network for the adaptive lasso weights as needed.
#
# Argument "data" must be a tibble with the following variables:
#   id
#   s_il
#   all elements of "vars"
#
# This function returns a list of 3 objects: 
#   weights,
#   estimated edges, 
#   edge_detection results (tp, tn, fp, fn, sse, etc.), 

EFL0_time_varying_graphical = function(data, 
                                       vars, 
                                       visits_of_interest, 
                                       visits_adaptive, 
                                       precision_edges_true = NULL,
                                       kernel = FALSE,
                                       nh = 9,
                                       h_power_low = -0.5,
                                       h_power_high = 1.5,
                                       adaptive = FALSE,
                                       adaptive_ridge = FALSE,
                                       L0 = FALSE,
                                       EFL0 = FALSE,
                                       beta_diff_max = 1,
                                       lambda2_seq_input = NULL,
                                       nfolds = 5,
                                       foldid_subj = NULL,
                                       seed = 123){
  
  set.seed(seed)
  
  # Specifications to add to output tibbles
  method = paste(ifelse(kernel == TRUE, "Kernel", "Visit-specific"), 
                 ifelse(adaptive == TRUE & adaptive_ridge == FALSE, " + Adaptive (OLS)", ""), 
                 ifelse(adaptive == TRUE & adaptive_ridge == TRUE,  " + Adaptive", ""),
                 ifelse(L0 == TRUE & EFL0 == FALSE, " + L0", ""), 
                 ifelse(L0 == TRUE & EFL0 == TRUE,  " + EFL0", ""),
                 sep = "")
  
  nvisits_of_interest = length(visits_of_interest)
  
  # Create foldid to divide SUBJECTS into k cross-validation folds
  if(is.null(foldid_subj)){
    
    ids = data %>%
      select(id) %>%
      distinct() %>%
      pull(id)
    
    if(length(ids) == 1){
      foldid_subj = data %>%
        select(id, s_il) %>%
        mutate(foldid = sample(rep(1:nfolds, length = nrow(data)), size = nrow(data), replace = FALSE))
    } else{
      foldid_subj = tibble(id = ids,
                           foldid = sample(rep(1:nfolds, length = length(ids)), size = length(ids), replace = FALSE))
    }
  } else{
    nfolds = nrow(foldid_subj %>% select(foldid) %>% distinct())
  }
  
  # Shells for analysis output
  visit_weights_all_h   = list()
  precision_edges_all_h = list()
  partial_edges_all_h   = list()
  edge_detection_all_h  = list()
  output_steps_all_h    = list()

  ################################################################################
  # Center each node on its visit-specific predicted mean using a loess regression
  
  loess_function = function(data){
    fit = loess(value ~ s_il, data = data)
    data %>%
      modelr::add_predictions(fit) %>%
      modelr::add_residuals(fit)
  }
  
  data_centered = data %>%
    pivot_longer(names_to = "key",
                 values_to = "value",
                 all_of(vars)) %>%
    select(id, s_il, key, value) %>%
    nest(data = c(id, value, s_il)) %>%
    mutate(data = map(data, loess_function)) %>%
    unnest(cols = data) %>%
    select(id, s_il, key, resid) %>%
    pivot_wider(names_from = key,
                values_from = resid)
  
  ################################################################################
  # Estimate the bandwidth sequence and loop through bandwidths to sum up cross-validation error
  # For visit-specific, need to make sure we're only using each observation for ONE time point's graph.
  #     Only look for a bandwidth as big as 1/2 the smallest difference between two visits.
  
  h_rule_of_thumb = density(x = data_centered$s_il, kernel = "gaussian", bw = "nrd0")$bw
  
  if(kernel){
    h_sequence = h_rule_of_thumb * exp(seq(from = h_power_low, to = h_power_high, length.out = nh))
  } else{
    diff_visits = NULL
    for(i in 2:(nvisits_of_interest)){
      diff_visits[i-1] = visits_of_interest[i] - visits_of_interest[i-1]
    }
    min_diff = min(diff_visits)
    h_sequence = min_diff / 2
  }
  
  print("h sequence:")
  print(h_sequence)
  
  cvm_bandwidth = rep(0, length(h_sequence))
  
  for(h in 1:length(h_sequence)){
    
    h_current = h_sequence[h]
    
    # Shells for analysis output
    visit_weights    = NULL
    precision_edges  = NULL
    partial_edges    = NULL
    edge_detection   = NULL
    adaptive_weights = NULL
    ridgeC           = NULL
    output_steps     = NULL
    
    ################################################################################
    # Loop through visits
    
    for(l in 1:nvisits_of_interest){
      
      s_l = visits_adaptive[l]
      
      orig_index = which(visits_of_interest == s_l)
      
      adaptive_index = ifelse(s_l < 0, which(visits_adaptive == visits_of_interest[[orig_index + 1]]),
                              ifelse(s_l > 0, which(visits_adaptive == visits_of_interest[[orig_index - 1]]),
                                     NA))
      
      ################################################################################
      # True edges
      
      if(!is.null(precision_edges_true)){
        precision_edges_true_current = precision_edges_true %>% 
          filter(s_il == s_l) %>% 
          select(row, col, edge) %>%
          arrange(row, col) 
      } else{
        precision_edges_true_current = NULL
      }
      
      ################################################################################
      # Construct smooth time-varying weights (either kernel or binary). 
      # For kernel, when s_l is negative (pre-diagnosis), only use pre-diagnosis observations.
      # Similarly, when s_l is positive (post-diagnosis), only use post-diagnosis observations.
      # For kernel, make the weights sum to the total number of observations (kernel * N / sum(kernel))
      # For binary, delete the rows with weight of 0
      
      if(kernel == TRUE){
        data_centered_visit = data_centered %>% 
          mutate(dist = (s_il - s_l)/h_current,
                 kernel_function = 1/h_current * (1/sqrt(2 * pi) * exp(-1/2 * dist^2))) %>%
          mutate(weight = kernel_function * n() / sum(kernel_function))
        
      } else{
        data_centered_visit = data_centered %>% 
          mutate(dist = (s_il - s_l)/h_current,
                 weight = as.numeric(-1 <= dist & dist < 1)) %>%
          filter(weight > 0)
      }
      
      ################################################################################
      # Run the EFL0_graphical() function for this time point
      
      if((adaptive == TRUE | EFL0 == TRUE) & l > 1){
        print(paste0("Time ", l, ": first case"))
        EFL0_results = EFL0_graphical(data_current = data_centered_visit,
                                      vars = vars,
                                      adaptive = adaptive,
                                      adaptive_current = adaptive,
                                      adaptive_current_wbeta = adaptive_weights[[adaptive_index]],
                                      adaptive_ridge = adaptive_ridge,
                                      L0 = L0,
                                      EFL0 = EFL0,
                                      EFL0_current = EFL0,
                                      EFL0_current_ridgeC = ridgeC[[adaptive_index]],
                                      beta_diff_max = beta_diff_max,
                                      lambda2_seq_input = lambda2_seq_input,
                                      seed = seed,
                                      nfolds = nfolds,
                                      foldid_subj = foldid_subj,
                                      precision_edges_true_current = precision_edges_true_current)
      } else{
        print(paste0("Time ", l, ": second case"))
        EFL0_results = EFL0_graphical(data_current = data_centered_visit,
                                      vars = vars,
                                      adaptive = adaptive,
                                      adaptive_current = FALSE,
                                      adaptive_ridge = adaptive_ridge,
                                      L0 = L0,
                                      EFL0 = EFL0,
                                      EFL0_current = FALSE,
                                      seed = seed,
                                      nfolds = nfolds,
                                      foldid_subj = foldid_subj,
                                      precision_edges_true_current = precision_edges_true_current)
      }
  
      visit_weights    = bind_rows(visit_weights,   data_centered_visit          %>% mutate(time = s_l, method = method) %>% select(time, s_il, weight, method))
      precision_edges  = bind_rows(precision_edges, EFL0_results$precision_edges %>% mutate(time = s_l, method = method))
      partial_edges    = bind_rows(partial_edges,   EFL0_results$partial_edges   %>% mutate(time = s_l, method = method))
      edge_detection   = bind_rows(edge_detection,  EFL0_results$edge_detection  %>% mutate(time = s_l, method = method))
      cvm_bandwidth[h] = cvm_bandwidth[h] +         as.numeric(EFL0_results$edge_detection$cvm)
      output_steps     = bind_rows(output_steps,    EFL0_results$output_steps    %>% mutate(time = s_l, method = method))
      
      ################################################################################
      # Store betas for Adaptive lasso at the next visit
      
      if(adaptive == TRUE){
        adaptive_weights[[l]] = 1/abs(EFL0_results$coefs_adaptive_weights)
      } else{
        adaptive_weights[[l]] = NULL
      }
      
      ################################################################################
      # Store ridgeC for EFL0 at the next visit
      
      if(EFL0 == TRUE){
        ridgeC[[l]] = EFL0_results$ridgeC_elastic_fuse
      } else{
        ridgeC[[l]] = NULL
      }
    }
    
    ################################################################################
    # Store results for this bandwidth in the lists of all results
    
    visit_weights_all_h[[h]]   = visit_weights
    precision_edges_all_h[[h]] = precision_edges
    partial_edges_all_h[[h]]   = partial_edges
    edge_detection_all_h[[h]]  = edge_detection
    output_steps_all_h[[h]]    = output_steps
  }
  
  ################################################################################
  # Select the bandwidth with smallest cvm and select the corresponding final results
  
  h_index = which(cvm_bandwidth == min(cvm_bandwidth))
  h_final = h_sequence[h_index]
  print(paste0("Final bandwidth h_index: ", h_index))
  
  visit_weights_final   = visit_weights_all_h[[h_index]]
  precision_edges_final = precision_edges_all_h[[h_index]]
  partial_edges_final   = partial_edges_all_h[[h_index]]
  edge_detection_final  = edge_detection_all_h[[h_index]] %>%
    mutate(h = h_final) %>%
    mutate(h_is_min_h = as.numeric(h == min(h_sequence)),
           h_is_max_h = as.numeric(h == max(h_sequence)))
  output_steps_final    = output_steps_all_h[[h_index]]
  
  return(list(visit_weights   = visit_weights_final,
              precision_edges = precision_edges_final,
              partial_edges   = partial_edges_final,
              edge_detection  = edge_detection_final,
              foldid          = foldid_subj,
              output_steps    = output_steps_final,
              h_selection     = tibble(h   = h_sequence,
                                       cvm = cvm_bandwidth)))
}