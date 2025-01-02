# Install and load required package
if (!requireNamespace("pso", quietly = TRUE)) install.packages("pso")
library(pso)

# LASSO Objective Function
lasso_objective <- function(beta, X, y, lambda) {
  n <- nrow(X)
  mse <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(mse + penalty)
}

# PSO for LASSO
pso_lasso <- function(X, y, lambda, swarm_size = 50, max_iter = 100) {
  p <- ncol(X)
  n <- nrow(X)
  
  # Define the function to minimize
  objective <- function(beta) lasso_objective(beta, X, y, lambda)
  
  # Run PSO
  result <- psoptim(
    par = rep(0, p),      # Initial beta
    fn = objective,       # Objective function
    lower = rep(-10, p),  # Lower bounds for beta
    upper = rep(10, p),   # Upper bounds for beta
    control = list(maxit = max_iter, s = swarm_size) # PSO parameters
  )
  
  # Return optimized beta
  return(result$par)
}

# Step 1: Define parameters
lambda <- best_lambda  # Use the optimal lambda from earlier
swarm_sizes <- c(10, 50, 100)
iterations <- c(50, 100, 200)

# Step 2: Evaluate impact of swarm size and iterations
results <- list()
for (s in swarm_sizes) {
  for (iter in iterations) {
    cat(sprintf("Swarm Size: %d, Iterations: %d\n", s, iter))
    beta_pso <- pso_lasso(X, y, lambda, swarm_size = s, max_iter = iter)
    results[[paste(s, iter, sep = "_")]] <- beta_pso
  }
}

# Compare results
for (key in names(results)) {
  cat(sprintf("Swarm Size & Iterations: %s\n", key))
  print(results[[key]])
}
