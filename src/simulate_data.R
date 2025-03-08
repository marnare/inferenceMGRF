simulate_outcomes <- function(N, k, n_outcomes = 2, rho = 0.5, seed = 1) {
  set.seed(seed)

  # Create single set of variances for all individuals
  variances <- numeric(n_outcomes)
  for(j in 1:n_outcomes) {
    variances[j] <- 0.1 - (j-1)*0.02  # Decreasing variance
  }

  # Create correlation matrix with constant correlation rho
  corr_matrix <- matrix(rho, nrow = n_outcomes, ncol = n_outcomes)
  diag(corr_matrix) <- 1

  # Create single covariance matrix for all individuals
  Sigma <- matrix(NA, nrow = n_outcomes, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    for(l in 1:n_outcomes) {
      Sigma[j,l] <- corr_matrix[j,l] * sqrt(variances[j] * variances[l])
    }
  }

  # Generate covariates (3 standard normal variables)
  X <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
  colnames(X) <- paste0("X", 1:3)

  # Generate true effects for each individual using the same covariance matrix
  ate_true <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(i in 1:N) {
    ate_true[i,] <- MASS::mvrnorm(1, mu = rep(5, n_outcomes), Sigma = Sigma)
  }

  # Treatment assignment
  P = rbinom(N, 1, 0.5)

  # Simulate outcomes
  Y <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    Y[,j] <- 0.5 + P*(ate_true[,j] + k) + rnorm(N, 0, 1)
  }

  # Create data frame with covariates
  data <- data.frame(P = P, X)

  # Add outcomes
  for(j in 1:n_outcomes) {
    data[[paste0("y", j)]] <- Y[,j]
  }

  # Add true effects
  for(j in 1:n_outcomes) {
    data[[paste0("ate", j, "_true")]] <- ate_true[,j] + k
    data[[paste0("ate", j, "_var")]] <- variances[j]  # Same variance for everyone
  }

  # Add covariances (same for everyone)
  for(j in 1:(n_outcomes-1)) {
    for(l in (j+1):n_outcomes) {
      data[[paste0("cov_", j, "_", l)]] <- Sigma[j,l]  # Same covariance for everyone
    }
  }

  return(data)
}



simulate_outcomes_treatments <- function(N, k, n_outcomes = 2, n_levels = 3, rho = 0.5, seed = 1) {
  set.seed(seed)
  
  # Create single set of variances for all individuals
  variances <- numeric(n_outcomes)
  for(j in 1:n_outcomes) {
    variances[j] <- 0.1 * exp(-0.1 * (j-1))
  }
  
  # Create correlation matrix for outcomes
  corr_matrix <- matrix(rho, nrow = n_outcomes, ncol = n_outcomes)
  diag(corr_matrix) <- 1
  
  # Create single covariance matrix for outcomes
  Sigma <- matrix(NA, nrow = n_outcomes, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    for(l in 1:n_outcomes) {
      Sigma[j,l] <- corr_matrix[j,l] * sqrt(variances[j] * variances[l])
    }
  }
  
  # Generate covariates (3 standard normal variables)
  X <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
  colnames(X) <- paste0("X", 1:3)
  
  # Generate SINGLE treatment with multiple levels
  P <- factor(sample(0:(n_levels-1), N, replace = TRUE))
  
  # Simulate outcomes with heterogeneous effects
  Y <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    # Base effect including covariate main effects
    Y[,j] <- 0.5 + 
      0.5 * X[,1] +    # Main effect of X1
      0.3 * X[,2] +    # Main effect of X2
      0.2 * X[,3] +    # Main effect of X3
      rnorm(N, 0, 1)
    
    # Add treatment effects with heterogeneity
    for(level in 1:n_levels) {
      level_indicator <- (as.numeric(as.character(P)) == (level - 1))
      
      # Treatment effect varies with covariates
      heterogeneous_effect <- level * (
        k +                    # Base treatment effect
          1.0 * X[,1] +         # Strong interaction with X1
          0.7 * X[,2] +         # Moderate interaction with X2
          0.5 * X[,3] +         # Weak interaction with X3
          0.3 * X[,1] * X[,2]   # Interaction between covariates
      )
      
      Y[,j] <- Y[,j] + level_indicator * heterogeneous_effect
    }
  }
  
  # Create data frame with only P, outcomes, and covariates
  data <- data.frame(
    P = P,
    X1 = X[,1],
    X2 = X[,2],
    X3 = X[,3]
  )
  
  # Add outcomes
  for(j in 1:n_outcomes) {
    data[[paste0("y", j)]] <- Y[,j]
  }
  
  return(data)
}





# And the multiple datasets simulation function:
simulate_multiple_datasets <- function(N = 500, k = 5, n_datasets = 500, 
                                       configs = list(
                                         list(outcomes = 2, levels = 2),
                                         list(outcomes = 10, levels = 3),
                                         list(outcomes = 30, levels = 3),
                                         list(outcomes = 50, levels = 3)
                                       ),
                                       rho = 0.5, output_path, name_prefix) {
  
  all_simulations <- list()
  
  for(config in configs) {
    n_outcomes <- config$outcomes
    n_levels <- config$levels
    
    cat(sprintf("\nSimulating %d datasets with %d outcomes and %d treatment levels\n", 
                n_datasets, n_outcomes, n_levels))
    
    pb <- txtProgressBar(min = 0, max = n_datasets, style = 3)
    config_sims <- list()
    
    for(i in 1:n_datasets) {
      config_sims[[i]] <- simulate_outcomes_treatments(
        N = N,
        k = k,
        n_outcomes = n_outcomes,
        n_levels = n_levels,
        rho = rho,
        seed = i
      )
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    config_name <- sprintf("outcomes_%d_levels_%d", n_outcomes, n_levels)
    all_simulations[[config_name]] <- config_sims
    
    saveRDS(all_simulations, 
            file = paste0(output_path, name_prefix, "simulated_datasets_temp.rds"))
  }
  
  saveRDS(all_simulations, 
          file = paste0(output_path, name_prefix, "simulated_datasets.rds"))
  
  return(all_simulations)
}






simulate_outcomes_personalized_varcov <- function(N, 
                                                  k, 
                                                  n_outcomes = 2, 
                                                  rho = 0.5, 
                                                  seed = 0) {
  
  set.seed(seed)
  
  # Generate covariates (3 standard normal variables)
  X <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
  colnames(X) <- paste0("X", 1:3)
  
  # Generate individual-specific variances for each outcome
  variances <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    base_var <- 1.0 - (j-1)*0.2  # Increased base variance
    variances[,j] <- rlnorm(N, meanlog = log(base_var), sdlog = 0.3)
  }
  
  # Create correlation matrix with constant correlation rho
  corr_matrix <- matrix(rho, nrow = n_outcomes, ncol = n_outcomes)
  diag(corr_matrix) <- 1
  
  # Generate true effects for each individual
  ate_true <- matrix(NA, nrow = N, ncol = n_outcomes)
  
  # Store covariances for each individual
  covariances <- array(NA, dim = c(N, n_outcomes, n_outcomes))
  
  for(i in 1:N) {
    # Create covariance matrix for individual i
    Sigma_i <- matrix(NA, nrow = n_outcomes, ncol = n_outcomes)
    for(j in 1:n_outcomes) {
      for(l in 1:n_outcomes) {
        Sigma_i[j,l] <- corr_matrix[j,l] * sqrt(variances[i,j] * variances[i,l])
        covariances[i,j,l] <- Sigma_i[j,l]  # Store the covariance
      }
    }
    
    # Generate correlated effects
    ate_i <- MASS::mvrnorm(1, mu = rep(5, n_outcomes), Sigma = Sigma_i)
    ate_true[i,] <- ate_i
  }
  
  # Treatment assignment
  P = rbinom(N, 1, 0.5)
  
  # Simulate outcomes with covariate effects
  Y <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    Y[,j] <- 0.5 + P*(ate_true[,j] + k) + 
      0.2*X[,1] + 0.3*X[,2] + 0.4*X[,3] +  # Add covariate effects
      rnorm(N, 0, 0.5)
  }
  
  # Create data frame with covariates
  data <- data.frame(P = P, X)
  
  # Add outcomes
  for(j in 1:n_outcomes) {
    data[[paste0("y", j)]] <- Y[,j]
  }
  
  # Add true effects and variances
  for(j in 1:n_outcomes) {
    data[[paste0("ate", j, "_true")]] <- ate_true[,j] + k
    data[[paste0("ate", j, "_var")]] <- variances[,j]
  }
  
  # Add covariances (using stored individual-specific covariances)
  for(j in 1:(n_outcomes-1)) {
    for(l in (j+1):n_outcomes) {
      data[[paste0("cov_", j, "_", l)]] <- covariances[,j,l]
    }
  }
  
  return(data)
}


# Define the function for creating a styled kable table
create_styled_table <- function(data, digits = 4, table_attr = "class='small-table'") {
  data %>%
    kable(digits = digits, table.attr = table_attr) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    scroll_box(width = "100%", height = "500px")  # You can adjust the height as needed
}

