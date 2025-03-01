simulate_outcomes <- function(N, 
                              n_outcomes, 
                              k = 0, 
                              seed = 0, 
                              rho = 0.1) {
  set.seed(seed)
  
  # Create correlation matrix with high correlation (rho)
  # Create matrix of size (n_outcomes + 1) to include X1
  dim <- n_outcomes + 1
  R <- matrix(rho, nrow = dim, ncol = dim)
  diag(R) <- 1
  
  # Define standard deviations: first for X1, then decreasing for outcomes
  sds <- c(2, seq(1, 0.5, length.out = n_outcomes))
  
  # Convert correlation to covariance matrix
  D <- diag(sds)
  Sigma <- D %*% R %*% D
  
  # Generate correlated normal variables
  X <- MASS::mvrnorm(N, mu = rep(0, dim), Sigma = Sigma)
  
  # Treatment assignment
  P <- rbinom(N, 1, 0.5)
  
  # Generate outcomes and true effects
  outcomes <- list()
  true_effects <- list()
  
  # Create outcomes and store true effects
  for(i in 1:n_outcomes) {
    # True effect is from the corresponding column in X (after X1)
    ate_true <- X[, i + 1]
    true_effects[[paste0("ate", i, "_true")]] <- ate_true + k
    
    # Generate outcome
    y <- 0.5 + P * (ate_true + k) + rnorm(N, 0, 0.5)
    outcomes[[paste0("y", i)]] <- y
  }
  
  # Combine all data
  data <- data.frame(
    do.call(cbind, outcomes),  # outcomes
    P = P,                     # treatment
    X = X[, 1],               # X1
    do.call(cbind, true_effects)  # true effects
  )
  
  return(data)
}






simulate_outcomes_personalized_varcov <- function(N, k, n_outcomes = 2, rho = 0.5) {
  set.seed(0)
  
  # Generate individual-specific variances for each outcome
  variances <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    base_var <- 0.1 - (j-1)*0.02  # Decreasing base variance
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
  
  # Simulate outcomes
  Y <- matrix(NA, nrow = N, ncol = n_outcomes)
  for(j in 1:n_outcomes) {
    Y[,j] <- 0.5 + P*(ate_true[,j] + k) + rnorm(N, 0, 0.5)
  }
  
  # Create data frame
  data <- data.frame(P = P)
  
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

