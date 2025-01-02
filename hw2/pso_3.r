# LASSO objective function
lasso_objective <- function(beta, X, y, lambda) {
  n <- nrow(X)
  mse <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(mse + penalty)
}

# PSO Algorithm
pso_lasso_custom <- function(X, y, lambda, swarm_size = 50, max_iter = 100, 
                              w = 0.5, c1 = 2, c2 = 2, lower_bound = -10, upper_bound = 10) {
  p <- ncol(X)
  # Initialize particles: positions and velocities
  positions <- matrix(runif(swarm_size * p, lower_bound, upper_bound), nrow = swarm_size, ncol = p)
  velocities <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  
  # Personal best positions and global best position
  personal_best <- positions
  personal_best_scores <- apply(personal_best, 1, function(beta) lasso_objective(beta, X, y, lambda))
  global_best <- personal_best[which.min(personal_best_scores), ]
  global_best_score <- min(personal_best_scores)
  
  # PSO iterations
  for (iter in 1:max_iter) {
    for (i in 1:swarm_size) {
      # Update velocity
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- w * velocities[i, ] +
        c1 * r1 * (personal_best[i, ] - positions[i, ]) +
        c2 * r2 * (global_best - positions[i, ])
      
      # Update position
      positions[i, ] <- positions[i, ] + velocities[i, ]
      
      # Enforce bounds
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)
      
      # Evaluate new position
      score <- lasso_objective(positions[i, ], X, y, lambda)
      
      # Update personal best
      if (score < personal_best_scores[i]) {
        personal_best[i, ] <- positions[i, ]
        personal_best_scores[i] <- score
      }
      
      # Update global best
      if (score < global_best_score) {
        global_best <- positions[i, ]
        global_best_score <- score
      }
    }
    
    # Optional: Print iteration progress
    cat(sprintf("Iteration %d, Best Score: %f\n", iter, global_best_score))
  }
  
  return(list(beta = global_best, score = global_best_score))
}

# Test PSO-LASSO
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)

lambda <- 0.1
result <- pso_lasso_custom(X, y, lambda, swarm_size = 50, max_iter = 1000)

# Print results
cat("Final LASSO Coefficients (PSO):\n")
print(result$beta)
