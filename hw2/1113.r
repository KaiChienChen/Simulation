coordinate_descent_lasso <- function(X, y, lambda, tol = 1e-6, max_iter = 1000) {
  n <- nrow(X) # 資料觀測值的數量
  p <- ncol(X) # 自變數的數量
  beta <- rep(0, p)  # 初始化 beta 為零向量，長度為 p
  beta_old <- beta # 儲存上一輪的 beta，用於收斂判斷
  iter <- 0
  converged <- FALSE # 初始化收斂標誌
  soft_threshold <- function(z, lambda) {
    sign(z) * pmax(0, abs(z) - lambda)
  }
  while (!converged && iter < max_iter) {
    for (j in 1:p) {
      # Compute the residual without the effect of the current beta_j
      r_j <- y - X %*% beta + X[, j] * beta[j]
      # Update beta_j
      rho <- sum(X[, j] * r_j)
      beta[j] <- soft_threshold(rho, lambda) / sum(X[, j]^2)
    }
    # Check convergence (if change in beta is small enough)
    if (sum(abs(beta - beta_old)) < tol) {
      converged <- TRUE
    }
    beta_old <- beta
    iter <- iter + 1
  }
  return(beta)
}
######## test
# j <- 1
# r_1 <- y - X %*% beta + X[, j] * beta[j]

########
# 使用範例：
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
head(X)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
lambda <- 0.1

# 跑LASSO
beta_est <- coordinate_descent_lasso(X, y , lambda)
print(beta_est)
