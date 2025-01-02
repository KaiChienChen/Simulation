# Required Libraries
library(glmnet)
# Step 1: Generate or load data
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
# Step 2: Standardize the predictors
X <- scale(X)

# Step 3: Define soft-thresholding operator
soft_threshold <- function(z, lambda) {
  sign(z) * pmax(0, abs(z) - lambda)
}

# Step 4: Coordinate Descent Function
coordinate_descent <- function(X, y, lambda, max_iter = 1000, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      beta[j] <- soft_threshold(sum(r_j * X[, j]) / n, lambda)
    }
    if (sqrt(sum((beta - beta_old)^2)) < tol) break
  }
  return(beta)
}

# Step 5: Cross-Validation for Lambda Selection
cv_lasso <- function(X, y, lambda_seq, folds = 5) {
  n <- nrow(X)
  fold_ids <- sample(rep(1:folds, length.out = n))
  cv_errors <- numeric(length(lambda_seq))
  
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    fold_errors <- numeric(folds)
    
    for (fold in 1:folds) {
      train_idx <- which(fold_ids != fold)
      test_idx <- which(fold_ids == fold)
      X_train <- X[train_idx, ]
      y_train <- y[train_idx]
      X_test <- X[test_idx, ]
      y_test <- y[test_idx]
      
      beta <- coordinate_descent(X_train, y_train, lambda)
      predictions <- X_test %*% beta
      fold_errors[fold] <- mean((y_test - predictions)^2)
    }
    cv_errors[i] <- mean(fold_errors)
  }
  
  best_lambda <- lambda_seq[which.min(cv_errors)]
  return(list(best_lambda = best_lambda, cv_errors = cv_errors))
}

# Step 6: Run Cross-Validation
lambda_seq <- seq(0.01, 5, length.out = 100)
cv_results <- cv_lasso(X, y, lambda_seq)

# Step 7: Compute Final Estimates with Best Lambda
best_lambda <- cv_results$best_lambda
final_beta <- coordinate_descent(X, y, best_lambda)

# Output Results
print(paste("Optimal Lambda:", best_lambda))
print("Final Coefficients:")
print(final_beta)


# Plot CV errors vs. Lambda
par(mfrow = c(1, 1))
plot(
  lambda_seq, cv_results$cv_errors, type = "l", col = "blue", lwd = 2,
  xlab = expression(lambda), ylab = "Cross-Validation Error",
  main = "Cross-Validation for LASSO"
)
abline(v = best_lambda, col = "red", lty = 2)
legend("topright", legend = c("CV Error", "Best Lambda"), col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# Compare true coefficients and estimated coefficients
par(mfrow = c(1, 2))

# True coefficients
plot(
  beta_true, type = "h", col = "darkgreen", lwd = 2,
  xlab = "Index", ylab = "Coefficient Value",
  main = "True Coefficients"
)

# Estimated coefficients
plot(
  final_beta, type = "h", col = "blue", lwd = 2,
  xlab = "Index", ylab = "Coefficient Value",
  main = "Estimated Coefficients"
)
