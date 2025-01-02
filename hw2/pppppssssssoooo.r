# R implementation of PSO for f(x, y) = x^2 + y^2
set.seed(123)

# Objective function
objective <- function(x, y) x^2 + y^2

# Parameters
N <- 20  # Number of particles
d <- 2   # Dimensions
iterations <- 100
omega <- 0.5
c1 <- 2
c2 <- 2

# Initialize particles
particles <- matrix(runif(N * d, -10, 10), ncol = d)
velocities <- matrix(runif(N * d, -1, 1), ncol = d)
personal_best <- particles
personal_best_value <- apply(personal_best, 1, function(p) objective(p[1], p[2]))
global_best <- personal_best[which.min(personal_best_value), ]
global_best_value <- min(personal_best_value)

# PSO loop
for (t in 1:iterations) {
  for (i in 1:N) {
    # Update velocity
    r1 <- runif(d)
    r2 <- runif(d)
    velocities[i, ] <- omega * velocities[i, ] +
      c1 * r1 * (personal_best[i, ] - particles[i, ]) +
      c2 * r2 * (global_best - particles[i, ])
    
    # Update position
    particles[i, ] <- particles[i, ] + velocities[i, ]
    
    # Evaluate
    value <- objective(particles[i, 1], particles[i, 2])
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
  
  cat("Iteration:", t, "Global Best Value:", global_best_value, "\n")
}

cat("Optimal solution found at:", global_best, "with value:", global_best_value, "\n")



