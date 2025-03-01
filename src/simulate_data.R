simulate_two_outcomes <- function(N, cov_dim, k){
  set.seed(0)
  
  # First create correlation matrix
  rho <- 0.9  # high correlation
  R <- matrix(c(1, rho, rho,
                rho, 1, rho,
                rho, rho, 1), nrow = 3)
  
  # Define standard deviations for each variable
  sds <- c(2, 1, 0.5)  # different scales for each variable
  
  # Convert correlation to covariance matrix
  # Sigma = D * R * D where D is diagonal matrix of standard deviations
  D <- diag(sds)
  Sigma <- D %*% R %*% D
  
  # Generate correlated normal variables
  X <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = Sigma)
  
  # Rest of the function remains the same...
  P = rbinom(N, 1, 0.5)
  
  ate1_true = X[, 2] 
  ate2_true = X[, 3]

  
  ## simulate outcome
  y1 <- 0.5 + P*(ate1_true + k) + rnorm(N, 0, 0.5)
  y2 <- 0.5 + P*(ate2_true + k) + rnorm(N, 0, 0.5)
  
  data = as.data.frame(cbind(y1, y2, P, X, 
                             ate1_true = (ate1_true + k), 
                             ate2_true = (ate2_true + k)))
  
  return(data)
}







simulate_two_outcomes_personalized_varcov <- function(N, cov_dim, k){
  set.seed(0)
  
  # Generate individual-specific variances and covariances
  # Using log-normal to ensure positive variances
  var1 <- rlnorm(N, meanlog = log(0.08), sdlog = 0.3)  # centered around 0.08
  var2 <- rlnorm(N, meanlog = log(0.06), sdlog = 0.3)  # centered around 0.06
  
  # Generate covariances that ensure positive definite matrices
  # Covariance needs to satisfy: |cov| <= sqrt(var1 * var2)
  max_covs <- sqrt(var1 * var2)
  # Generate correlation coefficients between 0.3 and 0.7
  corrs <- runif(N, 0.3, 0.7)
  covs <- corrs * max_covs
  
  # Generate true effects for each individual using their specific covariance matrix
  ate_true <- matrix(NA, nrow = N, ncol = 2)
  
  for(i in 1:N) {
    Sigma_i <- matrix(c(var1[i], covs[i],
                        covs[i], var2[i]), nrow = 2)
    # Generate correlated effects
    ate_i <- MASS::mvrnorm(1, mu = c(5, 5), Sigma = Sigma_i)
    ate_true[i,] <- ate_i
  }
  
  # Treatment assignment
  P = rbinom(N, 1, 0.5)
  
  # True effects
  ate1_true = ate_true[,1]
  ate2_true = ate_true[,2]
  
  ## simulate outcomes
  y1 <- 0.5 + P*(ate1_true + k) + rnorm(N, 0, 0.5)
  y2 <- 0.5 + P*(ate2_true + k) + rnorm(N, 0, 0.5)
  
  # Store individual variances and covariances
  data = data.frame(
    y1 = y1,
    y2 = y2,
    P = P,
    ate1_true = ate1_true + k,
    ate2_true = ate2_true + k,
    ate1_var = var1,
    ate2_var = var2,
    cov_1_2 = covs
  )
  
  return(data)
}










# Define the function for creating a styled kable table
create_styled_table <- function(data, digits = 4, table_attr = "class='small-table'") {
  data %>%
    kable(digits = digits, table.attr = table_attr) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    scroll_box(width = "100%", height = "500px")  # You can adjust the height as needed
}

