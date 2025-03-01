compute_bootstrap_stats <- function(X, 
                                    Y, 
                                    W, 
                                    num_iterations = 1000, 
                                    min.node.size = 50) {
  
  # Get number of outcomes from Y
  num_outcomes <- ncol(Y)
  
  # Initialize list to store predictions for each outcome
  predictions <- vector("list", num_outcomes)
  for(i in 1:num_outcomes) {
    predictions[[i]] <- list()
  }
  
  ## Bootstrapped predictions with a single tree
  for (i in 1:num_iterations) {
    model <- multi_arm_causal_forest(X, Y, W, 
                                     num.trees = 1, 
                                     seed = i, 
                                     min.node.size = min.node.size)
    
    # Store predictions for each outcome
    for(j in 1:num_outcomes) {
      predictions[[j]][[i]] <- model$predictions[, j]
    }
  }
  
  # Initialize results lists
  means <- vector("list", num_outcomes)
  vars <- vector("list", num_outcomes)
  
  # Compute means for each outcome
  for(i in 1:num_outcomes) {
    means[[i]] <- colMeans(do.call(rbind, predictions[[i]]), na.rm = TRUE)
  }
  
  # Initialize result dataframe
  result <- data.frame(matrix(nrow = nrow(X), ncol = 0))
  
  # Add means and compute variances
  for(i in 1:num_outcomes) {
    # Add means
    result[paste0("predictions", i, "_mean")] <- means[[i]]
    
    # Compute and add variances
    temp_var <- lapply(predictions[[i]], function(x) (x - means[[i]])^2)
    result[paste0("predictions", i, "_var")] <- colMeans(do.call(rbind, temp_var), na.rm = TRUE)
  }
  
  # Compute and add covariances
  for(i in 1:(num_outcomes-1)) {
    for(j in (i+1):num_outcomes) {
      temp_cov <- mapply(function(x, y) (x - means[[i]]) * (y - means[[j]]), 
                         predictions[[i]], predictions[[j]], SIMPLIFY = FALSE)
      result[paste0("cov_", i, "_", j)] <- colMeans(do.call(rbind, temp_cov), na.rm = TRUE)
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
    proportion_type1_risk = vol_only_nocov / mean(in_nocov),
    proportion_type2_risk = vol_only_cov / mean(in_cov),
    
    # Diagnostics
    dimension = d,
    n_samples = n_samples,
    mc_error = 1/sqrt(n_samples)
  )
  
  return(results)
}


