# LASSO objective function
lasso_objective <- function(beta, X, y, lambda) {
  n <- nrow(X)
  mse <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(mse + penalty)
}

# Custom PSO Implementation for LASSO
pso_lasso <- function(X, y, lambda, swarm_size = 100, max_iter = 1000, 
                      w = 0.5, c1 = 2, c2 = 2, lower_bound = -5, upper_bound = 5) {
  p <- ncol(X)
  
  # Initialize particles: positions and velocities
  positions <- matrix(runif(swarm_size * p, lower_bound, upper_bound), nrow = swarm_size, ncol = p)
  velocities <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  
  # Initialize personal bests and global best
  personal_best_positions <- positions
  personal_best_scores <- apply(personal_best_positions, 1, function(beta) lasso_objective(beta, X, y, lambda))
  global_best_position <- personal_best_positions[which.min(personal_best_scores), ]
  global_best_score <- min(personal_best_scores)
  fitness_values <- numeric(max_iter)
  # PSO iterations
  for (iter in 1:max_iter) {
    for (i in 1:swarm_size) {
      # Update velocity
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- w * velocities[i, ] +
        c1 * r1 * (personal_best_positions[i, ] - positions[i, ]) +
        c2 * r2 * (global_best_position - positions[i, ])
      
      # Update position
      positions[i, ] <- positions[i, ] + velocities[i, ]
      
      # Enforce bounds
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)
      
      # Evaluate fitness
      score <- lasso_objective(positions[i, ], X, y, lambda)
      
      # Update personal best
      if (score < personal_best_scores[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_scores[i] <- score
      }
      
      # Update global best
      if (score < global_best_score) {
        global_best_position <- positions[i, ]
        global_best_score <- score
      }
    }
    
    # Print progress
    cat(sprintf("Iteration %d, Best Score: %f\n", iter, global_best_score))
    fitness_values[iter] <- global_best_score
  }
  
  return(list(beta = global_best_position, score = global_best_score))
}

# Test PSO-LASSO Implementation
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)

lambda <- 0.5644444
result <- pso_lasso(X, y, lambda, swarm_size = 100, max_iter = 1000)

# Output results
cat("Optimized Coefficients (PSO-LASSO):\n")
print(result$beta)

# 保存每次迭代的最佳適應度
max_iter <- 1000
fitness_values <- numeric(max_iter)
for (iter in 1:max_iter) {
  # (在主迴圈內記錄 global_best_score)
  fitness_values[iter] <- global_best_score
}

# 繪製收斂曲線
plot(fitness_values, type = "l", col = "blue", lwd = 2, 
     xlab = "Iteration", ylab = "Best Fitness Value",
     main = "PSO Convergence")

############################################################3

# LASSO objective function
lasso_objective <- function(beta, X, y, lambda) {
  n <- nrow(X)
  mse <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(mse + penalty)
}

# Custom PSO Implementation with Particle Reset Mechanism
pso_lasso_with_reset <- function(X, y, lambda, swarm_size = 100, max_iter = 1000, 
                                 w = 0.9, c1 = 1, c2 = 1, 
                                 lower_bound = -5, upper_bound = 5, 
                                 reset_threshold = 10) {
  p <- ncol(X)
  
  # Initialize particles: positions, velocities, and reset counters
  positions <- matrix(runif(swarm_size * p, lower_bound, upper_bound), nrow = swarm_size, ncol = p)
  velocities <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  reset_counters <- rep(0, swarm_size)  # Track how long a particle is at the boundary
  # Initialize personal bests and global best
  personal_best_positions <- positions
  personal_best_scores <- apply(personal_best_positions, 1, function(beta) lasso_objective(beta, X, y, lambda))
  global_best_position <- personal_best_positions[which.min(personal_best_scores), ]
  global_best_score <- min(personal_best_scores)
  # PSO iterations
  for (iter in 1:max_iter) {
    w <- 0.9 - iter / max_iter * (0.9 - 0.4)
    for (i in 1:swarm_size) {
      # Update velocity
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- w * velocities[i, ] +
        c1 * r1 * (personal_best_positions[i, ] - positions[i, ]) +
        c2 * r2 * (global_best_position - positions[i, ])
      
      # Update position
      positions[i, ] <- positions[i, ] + velocities[i, ]
      
      # Enforce bounds
      at_boundary <- (positions[i, ] <= lower_bound | positions[i, ] >= upper_bound)
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)
      
      # Increment reset counter for particles stuck at boundary
      if (any(at_boundary)) {
        reset_counters[i] <- reset_counters[i] + 1
      } else {
        reset_counters[i] <- 0  # Reset counter if particle moves off boundary
      }
      
      # Reset particle if it stays at boundary for too long
      if (reset_counters[i] >= reset_threshold) {
        positions[i, ] <- runif(p, lower_bound, upper_bound)  # Randomize position
        velocities[i, ] <- runif(p, -1, 1)  # Randomize velocity
        reset_counters[i] <- 0  # Reset the counter
      }
      
      # Evaluate fitness
      score <- lasso_objective(positions[i, ], X, y, lambda)
      
      # Update personal best
      if (score < personal_best_scores[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_scores[i] <- score
      }
      
      # Update global best
      if (score < global_best_score) {
        global_best_position <- positions[i, ]
        global_best_score <- score
      }
    }
    
    # Print progress
    cat(sprintf("Iteration %d, Best Score: %f\n", iter, global_best_score))
  }
  
  return(list(beta = global_best_position, score = global_best_score))
}

# Example: Testing the improved PSO
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
lambda <- 0.5644444
result <- pso_lasso_with_reset(X, y, lambda, swarm_size = 200, max_iter = 10000)
# Output results
cat("Optimized Coefficients (PSO-LASSO with Reset):\n")
print(result$beta)


###############################################################
pso_lasso_improved <- function(X, y, lambda, swarm_size = 100, max_iter = 15000, 
                               w_start = 0.9, w_end = 0.4, 
                               c1 = 1, c2 = 1, 
                               lower_bound = -5, upper_bound = 5) {
  p <- ncol(X)
  
  # 初始化粒子
  positions <- matrix(runif(swarm_size * p, lower_bound, upper_bound), nrow = swarm_size, ncol = p)
  velocities <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  reset_counters <- rep(0, swarm_size)  # 重置計數器

  # 初始化個人最佳和全局最佳
  personal_best_positions <- positions
  personal_best_scores <- apply(personal_best_positions, 1, function(beta) lasso_objective(beta, X, y, lambda))
  global_best_position <- personal_best_positions[which.min(personal_best_scores), ]
  global_best_score <- min(personal_best_scores)
  
  # 主迴圈
  for (iter in 1:max_iter) {
    # 動態更新慣性權重
    w <- w_start - iter / max_iter * (w_start - w_end)
    
    for (i in 1:swarm_size) {
      # 更新速度
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- w * velocities[i, ] +
        c1 * r1 * (personal_best_positions[i, ] - positions[i, ]) +
        c2 * r2 * (global_best_position - positions[i, ])
      
      # 更新位置
      positions[i, ] <- positions[i, ] + velocities[i, ]
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)  # 確保邊界
      
      # 邊界檢查
      at_boundary <- (positions[i, ] == lower_bound | positions[i, ] == upper_bound)
      if (any(at_boundary)) reset_counters[i] <- reset_counters[i] + 1 else reset_counters[i] <- 0

      # 動態重置
      reset_threshold <- round(5 + iter / 1000)  # 動態閥值
      if (reset_counters[i] >= reset_threshold) {
        positions[i, ] <- global_best_position + runif(p, -0.1, 0.1)  # 附近隨機初始化
        velocities[i, ] <- runif(p, -0.5, 0.5)
        reset_counters[i] <- 0
      }

      # 適應度計算
      score <- lasso_objective(positions[i, ], X, y, lambda)
      
      # 更新個人最佳
      if (score < personal_best_scores[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_scores[i] <- score
      }
      
      # 更新全局最佳
      if (score < global_best_score) {
        global_best_position <- positions[i, ]
        global_best_score <- score
      }
    }
    
    # 可選：每 1000 次迭代打印一次進度
    if (iter %% 1000 == 0) {
      cat(sprintf("Iteration %d, Best Score: %f\n", iter, global_best_score))
    }
  }
  
  return(list(beta = global_best_position, score = global_best_score))
}
# Example: Testing the improved PSO
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
lambda <- 0.5644444
result <- pso_lasso_improved(X, y, lambda, swarm_size = 200, max_iter = 10000)
# Output results
cat("Optimized Coefficients (PSO-LASSO with Reset):\n")
print(result$beta)


#########################################################

pso_lasso_visualize <- function(X, y, lambda, swarm_size = 100, max_iter = 1000,
                                w_start = 0.9, w_end = 0.4,
                                c1 = 1, c2 = 1, lower_bound = -5, upper_bound = 5) {
  p <- ncol(X)
  positions <- matrix(runif(swarm_size * p, lower_bound, upper_bound), nrow = swarm_size, ncol = p)
  velocities <- matrix(runif(swarm_size * p, -1, 1), nrow = swarm_size, ncol = p)
  
  personal_best_positions <- positions
  personal_best_scores <- apply(personal_best_positions, 1, function(beta) lasso_objective(beta, X, y, lambda))
  global_best_position <- personal_best_positions[which.min(personal_best_scores), ]
  global_best_score <- min(personal_best_scores)
  
  # 收集數據
  best_scores <- numeric(max_iter)  # 每次迭代的全局最佳分數
  position_snapshots <- list()     # 保存粒子位置
  
  for (iter in 1:max_iter) {
    w <- w_start - iter / max_iter * (w_start - w_end)
    for (i in 1:swarm_size) {
      r1 <- runif(p)
      r2 <- runif(p)
      velocities[i, ] <- w * velocities[i, ] +
        c1 * r1 * (personal_best_positions[i, ] - positions[i, ]) +
        c2 * r2 * (global_best_position - positions[i, ])
      positions[i, ] <- positions[i, ] + velocities[i, ]
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)
      
      score <- lasso_objective(positions[i, ], X, y, lambda)
      if (score < personal_best_scores[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_scores[i] <- score
      }
      if (score < global_best_score) {
        global_best_position <- positions[i, ]
        global_best_score <- score
      }
    }
    
    # 收集當前迭代的數據
    best_scores[iter] <- global_best_score
    position_snapshots[[iter]] <- positions
  }
  
  return(list(best_scores = best_scores, position_snapshots = position_snapshots))
}

# 目標函數
lasso_objective <- function(beta, X, y, lambda) {
  n <- nrow(X)
  rss <- sum((y - X %*% beta)^2) / (2 * n)
  penalty <- lambda * sum(abs(beta))
  return(rss + penalty)
}

# 執行 PSO 並記錄數據
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
lambda <- 0.5644444

result <- pso_lasso_visualize(X, y, lambda, swarm_size = 200, max_iter = 10000)

# 視覺化
library(ggplot2)

# 1. 全局最佳分數收斂曲線
best_scores <- result$best_scores
df_scores <- data.frame(iteration = 1:length(best_scores), best_score = best_scores)

ggplot(df_scores, aes(x = iteration, y = best_score)) +
  geom_line(color = "blue") +
  labs(title = "Global Best Score Over Iterations", x = "Iteration", y = "Best Score") +
  theme_minimal()

# 2. 粒子位置的降維與分佈（使用 PCA）
install.packages("gganimate")
library(gganimate)

positions <- result$position_snapshots
iterations <- length(positions)

# 選擇最後 2 個特徵以進行簡化
df_particles <- do.call(rbind, lapply(1:iterations, function(i) {
  data.frame(
    x = positions[[i]][, 1],
    y = positions[[i]][, 2],
    iteration = i
  )
}))

ggplot(df_particles, aes(x = x, y = y, color = as.factor(iteration))) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_d() +
  labs(title = "Particle Positions Over Iterations", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()
