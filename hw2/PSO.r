set.seed(123)

# Generate data
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)

# LASSO objective function
lasso_objective <- function(beta, X, y, lambda) {
  n <- length(y)
  mse <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(mse + penalty)
}

# PSO parameters
N <- 50         # Number of particles
iterations <- 100
lambda <- 1     # LASSO penalty term
omega <- 0.5    # Inertia weight
c1 <- 1.5       # Cognitive coefficient
c2 <- 1.5       # Social coefficient

# Initialize particles
particles <- matrix(runif(N * p, -1, 1), nrow = N, ncol = p)
velocities <- matrix(runif(N * p, -0.1, 0.1), nrow = N, ncol = p)
personal_best <- particles
personal_best_value <- apply(particles, 1, function(beta) lasso_objective(beta, X, y, lambda))
global_best <- personal_best[which.min(personal_best_value), ]
global_best_value <- min(personal_best_value)

# PSO loop
for (t in 1:iterations) {
  for (i in 1:N) {
    # Update velocity
    r1 <- runif(p)
    r2 <- runif(p)
    velocities[i, ] <- omega * velocities[i, ] +
      c1 * r1 * (personal_best[i, ] - particles[i, ]) +
      c2 * r2 * (global_best - particles[i, ])
    
    # Update position
    particles[i, ] <- particles[i, ] + velocities[i, ]
    
    # Evaluate the objective function
    value <- lasso_objective(particles[i, ], X, y, lambda)
    if (value < personal_best_value[i]) {
      personal_best[i, ] <- particles[i, ]
      personal_best_value[i] <- value
    }
  }
  
  # Update global best
  current_global_best <- which.min(personal_best_value)
  if (personal_best_value[current_global_best] < global_best_value) {
    global_best <- personal_best[current_global_best, ]
    global_best_value <- personal_best_value[current_global_best]
  }
  
  # Output progress
  cat("Iteration:", t, "Global Best Value:", global_best_value, "\n")
}

# Final LASSO estimates
cat("LASSO estimates found by PSO:", global_best, "\n")


