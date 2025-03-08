compute_bootstrap_predictions <- function(X, Y, W, num_iterations = 1000) {
  num_outcomes <- ncol(Y)
  n_samples <- nrow(X)
  
  # Initialize lists for predictions
  predictions <- vector("list", num_outcomes)
  for(k in 1:num_outcomes) {
    predictions[[k]] <- list()
  }
  
  # Run bootstrap iterations
  for (i in 1:num_iterations) {
    model <- multi_arm_causal_forest(X, Y, W, num.trees = 1, seed = i)
    for(k in 1:num_outcomes) {
      predictions[[k]][[i]] <- model$predictions[, k]
    }
  }
  
  # Calculate means for each outcome
  pred_means <- vector("list", num_outcomes)
  for(k in 1:num_outcomes) {
    pred_means[[k]] <- colMeans(do.call(rbind, predictions[[k]]), na.rm = TRUE)
  }
  
  # Initialize lists for variances and covariances
  temp_var_lists <- vector("list", num_outcomes)
  temp_cov_lists <- vector("list", choose(num_outcomes, 2))
  
  # Calculate deviations for variances and covariances
  for (i in 1:num_iterations) {
    # Variances
    for(k in 1:num_outcomes) {
      if(i == 1) temp_var_lists[[k]] <- list()
      temp_var_lists[[k]][[i]] <- (predictions[[k]][[i]] - pred_means[[k]])^2
    }
    
    # Covariances
    cov_idx <- 1
    for(k in 1:(num_outcomes-1)) {
      for(m in (k+1):num_outcomes) {
        if(i == 1) temp_cov_lists[[cov_idx]] <- list()
        temp_cov_lists[[cov_idx]][[i]] <- 
          (predictions[[k]][[i]] - pred_means[[k]]) * 
          (predictions[[m]][[i]] - pred_means[[m]])
        cov_idx <- cov_idx + 1
      }
    }
  }
  
  # Calculate final estimates
  result <- data.frame(matrix(nrow = n_samples, ncol = 0))
  
  # Add means and variances
  for(k in 1:num_outcomes) {
    result[[paste0("predictions", k, "_mean")]] <- pred_means[[k]]
    result[[paste0("predictions", k, "_var")]] <- 
      colMeans(do.call(rbind, temp_var_lists[[k]]), na.rm = TRUE)
  }
  
  # Add covariances
  cov_idx <- 1
  for(k in 1:(num_outcomes-1)) {
    for(m in (k+1):num_outcomes) {
      result[[paste0("cov_", k, "_", m)]] <- 
        colMeans(do.call(rbind, temp_cov_lists[[cov_idx]]), na.rm = TRUE)
      cov_idx <- cov_idx + 1
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





process_all_simulations <- function(all_simulations, 
                                    num_iterations = 1000, 
                                    p_name = "P", 
                                    X_names = paste0("X.", 1:3)) {
  # Set up parallel processing
  plan(multisession, workers = availableCores() - 1)
  
  # Get first simulation to detect outcome names
  first_sim <- all_simulations[[1]][[1]]
  y_names <- grep("^y\\d+$", names(first_sim), value = TRUE)
  cat(sprintf("Detected outcomes: %s\n", paste(y_names, collapse = ", ")))
  
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
      Y <- as.matrix(sim_data[, y_names])
      W <- as.factor(sim_data[, p_name])
      X <- as.matrix(sim_data[, X_names])
      
      # Compute bootstrap stats
      result <- compute_bootstrap_predictions(
        X, Y, W, 
        num_iterations = num_iterations
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




### This is for personalized treatment effects
# Function to calculate volumes for each observation
calculate_volumes_individual <- function(data) {
  

  # Store results
  all_volumes <- list()
  
  # Loop through each observation
  for(i in 1:nrow(data)) {
    # Create individual covariance matrices
    cov_matrix <- matrix(c(
      data$predictions1_var[i], data$cov_1_2[i],
      data$cov_1_2[i], data$predictions2_var[i]
    ), nrow = 2)
    
    nocov_matrix <- matrix(c(
      data$predictions1_var[i], 0,
      0, data$predictions2_var[i]
    ), nrow = 2)
    
    # Get mean vector (in this case, the predictions for this observation)
    mean_vector <- c(data$predictions1_mean[i], data$predictions2_mean[i])
    
    # Calculate volumes for this observation
    vol_i <- calculate_volumes_hd(
      cov_matrix = cov_matrix,
      nocov_matrix = nocov_matrix,
      mean_vector = mean_vector,
      alpha = 0.05
    )
    
    all_volumes[[i]] <- vol_i
  }
  
  # Summarize results
  volume_ratios <- sapply(all_volumes, function(x) x$volume_ratio)
  type1_risks <- sapply(all_volumes, function(x) x$proportion_type1_risk)
  type2_risks <- sapply(all_volumes, function(x) x$proportion_type2_risk)
  
  summary_stats <- list(
    volume_ratio_mean = mean(volume_ratios),
    volume_ratio_sd = sd(volume_ratios),
    type1_risk_mean = mean(type1_risks),
    type1_risk_sd = sd(type1_risks),
    type2_risk_mean = mean(type2_risks),
    type2_risk_sd = sd(type2_risks),
    all_individual_results = all_volumes
  )
  
  return(summary_stats)
}



process_all_volumes_personalized <- function(all_results, chunk_size = 10) {
  # Set up parallel processing but with fewer workers to reduce memory load
  n_cores <- min(parallel::detectCores() - 1, 4)  # Limit max cores
  plan(multisession, workers = n_cores)
  
  all_volumes <- list()
  
  # Process each N-rho pair sequentially
  for(pair_name in names(all_results)) {
    cat(sprintf("\nProcessing volumes for %s\n", pair_name))
    
    pair_results <- all_results[[pair_name]]
    n_sims <- length(pair_results)
    n_chunks <- ceiling(n_sims/chunk_size)
    
    # Initialize list for this pair
    pair_volumes <- list()
    
    # Process in chunks
    for(chunk in 1:n_chunks) {
      cat(sprintf("\nProcessing chunk %d of %d\n", chunk, n_chunks))
      
      # Get chunk indices
      start_idx <- (chunk-1)*chunk_size + 1
      end_idx <- min(chunk*chunk_size, n_sims)
      chunk_indices <- start_idx:end_idx
      
      # Process chunk in parallel
      chunk_results <- future_lapply(chunk_indices, function(i) {
        tryCatch({
          calculate_volumes_individual(pair_results[[i]])
        }, error = function(e) {
          cat(sprintf("\nError processing simulation %d: %s\n", i, e$message))
          return(NULL)
        })
      }, future.seed = TRUE)
      
      # Add chunk results to pair volumes
      pair_volumes <- c(pair_volumes, chunk_results)
      
      # Intermediate save
      all_volumes[[pair_name]] <- pair_volumes
      saveRDS(all_volumes, paste0(data_path, name_prefix, "all_volumes_temp.rds"))
    }
    
    # Remove NULL results
    pair_volumes <- Filter(Negate(is.null), pair_volumes)
    
    if(length(pair_volumes) > 0) {
      all_volumes[[pair_name]] <- pair_volumes
      # Save after each pair is complete
      saveRDS(all_volumes, paste0(data_path, name_prefix, "all_volumes.rds"))
    }
  }
  
  # Return to sequential processing
  plan(sequential)
  
  return(all_volumes)
}



# Optimized version of calculate_volumes_individual
calculate_volumes_individual_optimized <- function(data) {
  n_obs <- nrow(data)
  
  # Pre-allocate vectors for results
  volume_ratios <- numeric(n_obs)
  type1_risks <- numeric(n_obs)
  type2_risks <- numeric(n_obs)
  
  # Vectorize matrix creation
  cov_matrices <- array(0, dim = c(2, 2, n_obs))
  nocov_matrices <- array(0, dim = c(2, 2, n_obs))
  mean_vectors <- matrix(0, nrow = n_obs, ncol = 2)
  
  # Fill arrays efficiently
  cov_matrices[1,1,] <- data$predictions1_var
  cov_matrices[1,2,] <- data$cov_1_2
  cov_matrices[2,1,] <- data$cov_1_2
  cov_matrices[2,2,] <- data$predictions2_var
  
  nocov_matrices[1,1,] <- data$predictions1_var
  nocov_matrices[2,2,] <- data$predictions2_var
  
  mean_vectors[,1] <- data$predictions1_mean
  mean_vectors[,2] <- data$predictions2_mean
  
  # Process in chunks for memory efficiency
  chunk_size <- 1000
  n_chunks <- ceiling(n_obs/chunk_size)
  
  for(chunk in 1:n_chunks) {
    start_idx <- (chunk-1)*chunk_size + 1
    end_idx <- min(chunk*chunk_size, n_obs)
    chunk_indices <- start_idx:end_idx
    
    # Process chunk in parallel
    chunk_results <- future_lapply(chunk_indices, function(i) {
      calculate_volumes_hd(
        cov_matrix = cov_matrices[,,i],
        nocov_matrix = nocov_matrices[,,i],
        mean_vector = mean_vectors[i,],
        alpha = 0.05
      )
    }, future.seed = TRUE)
    
    # Extract results
    volume_ratios[chunk_indices] <- sapply(chunk_results, function(x) x$volume_ratio)
    type1_risks[chunk_indices] <- sapply(chunk_results, function(x) x$proportion_type1_risk)
    type2_risks[chunk_indices] <- sapply(chunk_results, function(x) x$proportion_type2_risk)
  }
  
  # Create summary statistics
  summary_stats <- list(
    volume_ratio_mean = mean(volume_ratios),
    volume_ratio_sd = sd(volume_ratios),
    type1_risk_mean = mean(type1_risks),
    type1_risk_sd = sd(type1_risks),
    type2_risk_mean = mean(type2_risks),
    type2_risk_sd = sd(type2_risks),
    all_individual_results = list(
      volume_ratios = volume_ratios,
      type1_risks = type1_risks,
      type2_risks = type2_risks
    )
  )
  
  return(summary_stats)
}




#### With many treatments

compute_bootstrap_predictions_multivalued <- function(X, Y, W, num_iterations = 1000) {
  num_outcomes <- ncol(Y)
  n_samples <- nrow(X)
  
  # Initialize lists for predictions
  predictions <- vector("list", num_outcomes)
  for(k in 1:num_outcomes) {
    predictions[[k]] <- list()
  }
  
  # Run bootstrap iterations
  for (i in 1:num_iterations) {
    model <- multi_arm_causal_forest(X, Y, W, num.trees = 1, seed = i)
    for(k in 1:num_outcomes) {
      predictions[[k]][[i]] <- model$predictions[, k]
    }
  }
  
  # Calculate means for each outcome
  pred_df <- data.frame(matrix(NA, nrow = n_samples, ncol = num_outcomes))
  colnames(pred_df) <- paste0("pred_y", 1:num_outcomes)
  
  for(k in 1:num_outcomes) {
    pred_df[,k] <- colMeans(do.call(rbind, predictions[[k]]), na.rm = TRUE)
  }
  
  # Calculate variance-covariance matrix for each observation
  vcov_matrices <- cov(pred_df)
  
  return(list(
    predictions = pred_df,
    vcov_matrices = vcov_matrices
  ))
}



process_all_simulations_multivalued <- function(all_simulations, 
                                                num_iterations = 1000, 
                                                p_name = "P", 
                                                X_names = paste0("X", 1:3)) {
  plan(multisession, workers = availableCores() - 1)
  
  all_results <- lapply(names(all_simulations), function(pair_name) {
    # Get y_names for this specific dataset
    first_sim_in_pair <- all_simulations[[pair_name]][[1]]
    y_names <- grep("^y\\d+$", names(first_sim_in_pair), value = TRUE)
    cat(sprintf("\nProcessing %s with outcomes: %s\n", 
                pair_name, paste(y_names, collapse = ", ")))
    
    n_sims <- length(all_simulations[[pair_name]])
    pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
    
    pair_results <- lapply(seq_len(n_sims), function(i) {
      sim_data <- all_simulations[[pair_name]][[i]]
      
      Y <- as.matrix(sim_data[, y_names])
      W <- as.factor(sim_data[, p_name])
      X <- as.matrix(sim_data[, X_names])
      
      result <- compute_bootstrap_predictions_multivalued(
        X, Y, W, 
        num_iterations = num_iterations
      )
      
      setTxtProgressBar(pb, i)
      return(result)
    })
    
    close(pb)
    
    all_results_temp <- setNames(list(pair_results), pair_name)
    saveRDS(all_results_temp, 
            paste0(data_path, name_prefix, "results_", pair_name, ".rds"))
    
    return(pair_results)
  })
  
  names(all_results) <- names(all_simulations)
  plan(sequential)
  
  saveRDS(all_results, paste0(data_path, name_prefix, "all_results.rds"))
  
  return(all_results)
}


process_volumes_multivalued <- function(predictions, print = FALSE) {
  # Get number of outcomes from the predictions data frame
  num_outcomes <- ncol(predictions$predictions)
  
  # Calculate means for all outcomes
  means <- colMeans(predictions$predictions)
  
  # Calculate actual covariance matrix between all predictions
  cov_matrix <- predictions$vcov_matrices
  
  # Create shape matrix (diagonal matrix with same variances)
  cov_matrix_shape <- matrix(0, nrow = num_outcomes, ncol = num_outcomes)
  diag(cov_matrix_shape) <- diag(cov_matrix)
  
  if(print) {
    cat("\nCovariance Matrix:\n")
    print(cov_matrix)
    cat("\nShape Matrix:\n")
    print(cov_matrix_shape)
    cat("\nMeans:\n")
    print(means)
  }
  
  # Calculate volumes
  volumes <- calculate_volumes_hd_multivalued(
    cov_matrix = cov_matrix,
    nocov_matrix = cov_matrix_shape,
    mean_vector = means
  )
  
  return(volumes)
}



process_all_volumes_multivalued <- function(all_results) {
  all_volumes <- list()
  
  for(pair_name in names(all_results)) {
    cat(sprintf("\nProcessing volumes for %s\n", pair_name))
    
    # Process each simulation in this pair
    pair_volumes <- lapply(all_results[[pair_name]], function(sim_result) {
      tryCatch({
        process_volumes_multivalued(sim_result)
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



calculate_volumes_hd_multivalued <- function(cov_matrix, nocov_matrix, mean_vector, alpha = 0.05, n_samples = 100000) {
  d <- length(mean_vector)
  chi_sq_crit <- qchisq(1 - alpha, df = d)
  
  # Monte Carlo estimation
  points_cov <- MASS::mvrnorm(n_samples, rep(0, d), cov_matrix)
  points_nocov <- MASS::mvrnorm(n_samples, rep(0, d), nocov_matrix)
  
  # Calculate Mahalanobis distances
  maha_cov_under_cov <- mahalanobis(points_cov, rep(0, d), cov_matrix)
  maha_cov_under_nocov <- mahalanobis(points_cov, rep(0, d), nocov_matrix)
  maha_nocov_under_cov <- mahalanobis(points_nocov, rep(0, d), cov_matrix)
  maha_nocov_under_nocov <- mahalanobis(points_nocov, rep(0, d), nocov_matrix)
  
  # Calculate empirical quantiles instead of using chi-square critical value
  quant_cov <- quantile(maha_cov_under_cov, 1 - alpha)
  quant_nocov <- quantile(maha_nocov_under_nocov, 1 - alpha)
  
  # Calculate risks using empirical distributions
  type1_risk <- mean(maha_cov_under_nocov > quant_nocov)
  type2_risk <- mean(maha_nocov_under_cov > quant_cov)
  
  results <- list(
    proportion_type1_risk = type1_risk,
    proportion_type2_risk = type2_risk,
    dimension = d,
    n_samples = n_samples,
    mc_error = 1/sqrt(n_samples),
    quantile_cov = quant_cov,
    quantile_nocov = quant_nocov
  )
  
  return(results)
}