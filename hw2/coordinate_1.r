coordinate_descent_lasso <- function(X, y, lambda, tol = 1e-6, max_iter = 1000) {
  n <- nrow(X) # 資料觀測值的數量
  p <- ncol(X) # 自變數的數量
  beta <- rep(0, p)  # 初始化 beta 為零向量
  beta_old <- beta # 儲存上一輪的 beta，用於收斂判斷
  iter <- 0
  converged <- FALSE # 初始化收斂標誌

  # 定義軟閾值函數
  soft_threshold <- function(z, lambda) {
    sign(z) * pmax(0, abs(z) - lambda)
  }

  while (!converged && iter < max_iter) {
    for (j in 1:p) {
      # 計算殘差不包括 beta_j 的影響
      r_j <- y - X %*% beta + X[, j] * beta[j]
      # 計算 rho
      rho <- sum(X[, j] * r_j)
      # 更新 beta_j
      beta[j] <- soft_threshold(rho, lambda) / sum(X[, j]^2)
    }

    # 判斷收斂條件
    if (sum(abs(beta - beta_old)) / sum(abs(beta_old) + tol) < tol) {
      converged <- TRUE
    }

    beta_old <- beta
    iter <- iter + 1
  }
  return(beta)
}

# 測試 LASSO
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
lambda <- 0.1

# 執行 coordinate descent LASSO
beta_est <- coordinate_descent_lasso(X, y, lambda)
print(beta_est)

#####################################################################

cv_lasso <- function(X, y, lambda_seq, K = 5, tol = 1e-6, max_iter = 1000) {
  n <- nrow(X)
  folds <- sample(rep(1:K, length.out = n)) # 隨機分成 K 個折
  mse <- numeric(length(lambda_seq)) # 儲存每個 lambda 的交叉驗證誤差

  # 對每個 lambda 進行交叉驗證
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    mse_fold <- numeric(K) # 儲存每個折的 MSE
    for (k in 1:K) {
      # 分割資料為訓練集和驗證集
      train_idx <- which(folds != k)
      test_idx <- which(folds == k)
      X_train <- X[train_idx, ]
      y_train <- y[train_idx]
      X_test <- X[test_idx, ]
      y_test <- y[test_idx]
      
      # 訓練 LASSO 模型
      beta <- coordinate_descent_lasso(X_train, y_train, lambda, tol, max_iter)
      
      # 預測和計算 MSE
      y_pred <- X_test %*% beta
      mse_fold[k] <- mean((y_test - y_pred)^2)
    }
    mse[i] <- mean(mse_fold) # 計算此 lambda 的平均 MSE
  }
  
  # 找到最佳 lambda
  best_lambda <- lambda_seq[which.min(mse)]
  
  return(list(best_lambda = best_lambda, mse = mse, lambda_seq = lambda_seq))
}

# 測試交叉驗證
set.seed(123)
lambda_seq <- 10^seq(0, 1, length.out = 10) # 在 0.01 到 10 的範圍內取對數序列
cv_result <- cv_lasso(X, y, lambda_seq)

# 顯示最佳 lambda 和對應的誤差
print(cv_result$best_lambda)

# 繪製 lambda 與 MSE 的圖
plot(log10(cv_result$lambda_seq), cv_result$mse, type = "b",
     xlab = "log10(Lambda)", ylab = "Mean Squared Error (MSE)",
     main = "Cross-Validation for LASSO")
abline(v = log10(cv_result$best_lambda), col = "red", lty = 2)
