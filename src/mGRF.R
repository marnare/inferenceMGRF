compute_bootstrap_stats <- function(X, Y, W, num_iterations = 1000, min.node.size = 50) {
  num_outcomes <- ncol(Y)
  
  # Use future_lapply with proper random seed handling
  predictions <- future_lapply(1:num_iterations, 
                               function(i) {
                                 # Create bootstrap sample
                                 boot_idx <- sample(1:nrow(X), replace = TRUE)
                                 X_boot <- X[boot_idx,]
                                 Y_boot <- Y[boot_idx,]
                                 W_boot <- W[boot_idx]
                                 
                                 # Fit model on bootstrap sample
                                 model <- multi_arm_causal_forest(X_boot, Y_boot, W_boot, 
                                                                  num.trees = 1, 
                                                                  seed = i, 
                                                                  min.node.size = min.node.size)
                                 
                                 # Predict on original data
                                 return(predict(model, X))
                               }, 
                               future.seed = TRUE)  # Add proper seed handling for parallel processing
  
  # Rest of the function remains the same
  pred_array <- array(unlist(predictions), 
                      dim = c(nrow(X), num_outcomes, num_iterations))
  
  # Initialize result dataframe
  result <- data.frame(matrix(nrow = nrow(X), ncol = 0))
  
  # Calculate means, variances, and covariances for each observation
  for(obs in 1:nrow(X)) {
    # Extract predictions for this observation
    obs_preds <- matrix(pred_array[obs,,], ncol = num_outcomes)
    
    # Calculate means
    means <- colMeans(obs_preds, na.rm = TRUE)
    
    # Calculate covariance matrix
    cov_matrix <- cov(obs_preds, use = "pairwise.complete.obs")
    
    # Store means
    for(i in 1:num_outcomes) {
      result[obs, paste0("predictions", i, "_mean")] <- means[i]
      result[obs, paste0("predictions", i, "_var")] <- cov_matrix[i,i]
    }
    
    # Store covariances
    for(i in 1:(num_outcomes-1)) {
      for(j in (i+1):num_outcomes) {
        result[obs, paste0("cov_", i, "_", j)] <- cov_matrix[i,j]
      }
    }
  }
  
  return(result)
}




calculate_volumes_hd <- function(cov_matrix, nocov_matrix, mean_vector, alpha = 0.05, n_samples = 100000) {
  d <- length(mean_vector)
  chi_sq_crit <- qchisq(1 - alpha, df = d)
  
  # Theoretical volumes using determinants
  # Volume of d-dimensional ellipsoid = (pi^(d/2)/gamma(d/2+1)) * sqrt(det(Sigma)) * radius^d
  volume_constant <- (pi^(d/2)/gamma(d/2 + 1)) * sqrt(chi_sq_crit)^d
  vol_cov_theoretical <- volume_constant * sqrt(det(cov_matrix))
  vol_nocov_theoretical <- volume_constant * sqrt(det(nocov_matrix))
  
  # Monte Carlo estimation
  # Generate random points from a large enough box containing both ellipsoids
  points <- MASS::mvrnorm(n_samples, rep(0, d), diag(d))
  
  
  # Check which points fall in each region using Mahalanobis distance
  in_cov <- mahalanobis(points, rep(0, d), cov_matrix) <= chi_sq_crit
  in_nocov <- mahalanobis(points, rep(0, d), nocov_matrix) <= chi_sq_crit
  
  # Calculate volumes and proportions
  vol_intersection <- mean(in_cov & in_nocov)
  vol_only_cov <- mean(in_cov & !in_nocov)
  vol_only_nocov <- mean(!in_cov & in_nocov)
  
  # Scale Monte Carlo proportions to actual volumes
  scaling_factor <- vol_cov_theoretical / mean(in_cov)
  
  results <- list(
    # Theoretical volumes
    volume_with_cov = vol_cov_theoretical,
    volume_without_cov = vol_nocov_theoretical,
    volume_ratio = vol_cov_theoretical / vol_nocov_theoretical,
    
    # Monte Carlo estimates
    volume_intersection = vol_intersection * scaling_factor,
    volume_only_cov = vol_only_cov * scaling_factor,
    volume_only_nocov = vol_only_nocov * scaling_factor,
    
    # Proportions
    proportion_intersection = vol_intersection / min(mean(in_cov), mean(in_nocov)),
    proportion_type2_risk = vol_only_nocov / mean(in_nocov),
    proportion_type1_risk = vol_only_cov / mean(in_cov),
    
    # Diagnostics
    dimension = d,
    n_samples = n_samples,
    mc_error = 1/sqrt(n_samples)
  )
  
  return(results)
}





## This is for average treatment effects
process_all_simulations <- function(all_simulations, num_iterations = 1000, min.node.size = 5) {
  # Set up parallel processing
  plan(multisession, workers = availableCores() - 1)
  
  # Process each N-rho combination using lapply
  all_results <- lapply(names(all_simulations), function(pair_name) {
    cat(sprintf("\nProcessing %s\n", pair_name))
    
    # Set up progress bar
    n_sims <- length(all_simulations[[pair_name]])
    pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
    
    # Process each simulation in this pair using lapply
    pair_results <- lapply(seq_len(n_sims), function(i) {
      sim_data <- all_simulations[[pair_name]][[i]]
      
      # Prepare data
      Y <- as.matrix(sim_data[, c(1, 2)])
      W <- as.factor(sim_data[, 3])
      X <- as.matrix(sim_data[, c(4:ncol(sim_data))])
      
      # Compute bootstrap stats
      result <- compute_bootstrap_stats(
        X, Y, W, 
        num_iterations = num_iterations,
        min.node.size = min.node.size
      )
      
      setTxtProgressBar(pb, i)
      return(result)
    })
    
    close(pb)
    
    # Save intermediate results
    all_results_temp <- setNames(list(pair_results), pair_name)
    saveRDS(all_results_temp, 
            paste0(data_path, name_prefix, "results_", pair_name, ".rds"))
    
    return(pair_results)
  })
  
  # Name the results list
  names(all_results) <- names(all_simulations)
  
  # Return to sequential processing
  plan(sequential)
  
  # Save final complete results
  saveRDS(all_results, paste0(data_path, name_prefix, "all_results.rds"))
  
  return(all_results)
}




## This is for average treatment effects
process_volumes <- function(predictions, print = FALSE) {
  # Calculate actual covariance between predictions
  pred_cov <- cov(predictions$predictions1_mean, predictions$predictions2_mean)
  
  # Create covariance matrices (population average treatment effects)
  cov_matrix <- matrix(c(
    var(predictions$predictions1_mean), 
    pred_cov,
    pred_cov, 
    var(predictions$predictions2_mean)
  ), nrow = 2, ncol = 2)
  
  # Create shape matrix (diagonal matrix with same variances)
  cov_matrix_shape <- matrix(c(
    var(predictions$predictions1_mean), 
    0,
    0, 
    var(predictions$predictions2_mean)
  ), nrow = 2, ncol = 2)
  
  if(print) {
    cat("\nCovariance Matrix:\n")
    print(cov_matrix)
    cat("\nShape Matrix:\n")
    print(cov_matrix_shape)
  }
  
  # Calculate volumes
  volumes <- calculate_volumes_hd(
    cov_matrix = cov_matrix,
    nocov_matrix = cov_matrix_shape,
    mean_vector = c(mean(predictions$predictions1_mean), 
                    mean(predictions$predictions2_mean))
  )
  
  return(volumes)
}



## This is for average treatment effects

# Modified process_all_volumes function
process_all_volumes <- function(all_results) {
  all_volumes <- list()
  
  for(pair_name in names(all_results)) {
    cat(sprintf("\nProcessing volumes for %s\n", pair_name))
    
    # Process each simulation in this pair
    pair_volumes <- lapply(all_results[[pair_name]], function(sim_result) {
      tryCatch({
        process_volumes(sim_result)
      }, error = function(e) {
        cat(sprintf("\nError processing simulation in %s:\n", pair_name))
        print(e)
        return(NULL)
      })
    })
    
    # Remove NULL results
    pair_volumes <- Filter(Negate(is.null), pair_volumes)
    
    if(length(pair_volumes) > 0) {
      all_volumes[[pair_name]] <- pair_volumes
      saveRDS(all_volumes, paste0(data_path, name_prefix, "all_volumes.rds"))
    }
  }
  
  return(all_volumes)
}




