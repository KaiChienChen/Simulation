# Load necessary libraries
# install.packages("glmnet")
library(glmnet)

# Define LASSO using coordinate descent
lasso_coord_descent <- function(X, Y, lambda, tol = 1e-6, max_iter = 1000) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)  # Initialize coefficients
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      # Compute partial residual
      r <- Y - X %*% beta + X[, j] * beta[j]
      # Update beta[j] with soft-thresholding
      beta[j] <- soft_threshold(sum(r * X[, j]) / n, lambda)
    }
    # Check for convergence
    if (sum(abs(beta - beta_old)) < tol) break
  }
  return(beta)
}

# Soft-thresholding function
soft_threshold <- function(z, lambda) {
  sign(z) * max(abs(z) - lambda, 0)
}

# Simulation setup
simulate_data <- function(n, p, beta_true, sigma = 1) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- X %*% beta_true + rnorm(n, sd = sigma)
  return(list(X = X, Y = Y))
}




# Define PSO for LASSO
pso_lasso <- function(X, Y, lambda, swarm_size = 30, max_iter = 100) {
  p <- ncol(X)
  # Initialize positions and velocities
  positions <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  velocities <- matrix(0, nrow = swarm_size, ncol = p)
  personal_best <- positions
  global_best <- positions[1, ]
  # Evaluate initial fitness
  fitness <- apply(positions, 1, function(beta) sum((Y - X %*% beta)^2) + lambda * sum(abs(beta)))
  personal_best_fitness <- fitness
  global_best_fitness <- min(fitness)
  for (iter in 1:max_iter) {
    # Update velocities and positions
    for (i in 1:swarm_size) {
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- 0.7 * velocities[i, ] +
        1.5 * r1 * (personal_best[i, ] - positions[i, ]) +
        1.5 * r2 * (global_best - positions[i, ])
      positions[i, ] <- positions[i, ] + velocities[i, ]
    }
    
    # Evaluate fitness and update personal/global bests
    fitness <- apply(positions, 1, function(beta) sum((Y - X %*% beta)^2) + lambda * sum(abs(beta)))
    for (i in 1:swarm_size) {
      if (fitness[i] < personal_best_fitness[i]) {
        personal_best[i, ] <- positions[i, ]
        personal_best_fitness[i] <- fitness[i]
      }
    }
    if (min(fitness) < global_best_fitness) {
      global_best <- positions[which.min(fitness), ]
      global_best_fitness <- min(fitness)
    }
  }
  return(global_best)
}






# Example usage
set.seed(123)
n <- 50
p <- 200
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
sim_data <- simulate_data(n, p, beta_true)
head(sim_data)
# Apply Coordinate Descent
lambda <- 0.1
beta_cd <- lasso_coord_descent(sim_data$X, sim_data$Y, lambda)
head(beta_cd)
# Apply PSO
beta_pso <- pso_lasso(sim_data$X, sim_data$Y, lambda)

# Evaluate performance
mse <- function(beta_est, beta_true) {
  mean((beta_est - beta_true)^2)
}
cat("MSE (CD):", mse(beta_cd, beta_true), "\n")
cat("MSE (PSO):", mse(beta_pso, beta_true), "\n")
