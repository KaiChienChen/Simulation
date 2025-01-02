# 定義協調下降函數
coordinate_descent <- function(X, y, lambda, alpha = 1, tol = 1e-4, max_iter = 2000) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p) # 初始化 beta 係數向量
  intercept <- mean(y) # 初始化截距項

  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      # 計算部分殘差
      residuals <- y - intercept - X %*% beta
      # 計算 j 處的梯度
      gradient <- -t(X[, j]) %*% residuals / n
      # 應用軟閾值
      beta[j] <- soft_threshold(gradient, lambda * alpha) / (1 + lambda * (1 - alpha))
    }
    # 檢查收斂性
    if (sum(abs(beta - beta_old)) < tol) {
      break
    }
  }
  return(list(intercept = intercept, beta = beta))
}

# 定義軟閾值函數
soft_threshold <- function(z, gamma) {
  return(sign(z) * max(0, abs(z) - gamma))
}

# 進行交叉驗證
cross_validation <- function(X, y, lambda_seq, alpha = 1, nfolds = 10) {
  n <- nrow(X)
  foldid <- sample(rep(1:nfolds, length.out = n)) # 建立隨機的折疊索引
  cv_errors <- numeric(length(lambda_seq))
  
  for (k in 1:length(lambda_seq)) {
    lambda <- lambda_seq[k]
    for (i in 1:nfolds) {
      # 將資料分為訓練集和測試集
      X_train <- X[foldid != i, ]
      y_train <- y[foldid != i]
      X_test <- X[foldid == i, ]
      y_test <- y[foldid == i]
      # 使用協調下降法擬合模型
      model <- coordinate_descent(X_train, y_train, lambda, alpha)
      # 計算測試集上的預測誤差
      y_pred <- model$intercept + X_test %*% model$beta
      cv_errors[k] <- cv_errors[k] + mean((y_test - y_pred)^2)
    }
    cv_errors[k] <- cv_errors[k] / nfolds
  }
  return(cv_errors)
}

################################
# 範例使用
# 生成一些範例資料
set.seed(123)
n <- 50
p <- 200
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3.5, 4, -2.8, 3.2, rep(0, p - 5))
y <- X %*% beta_true + rnorm(n)
# 定義 lambda 序列
lambda_seq <- exp(seq(0, 1, length.out = 10))

# 執行交叉驗證
cv_errors <- cross_validation(X, y, lambda_seq)

# 找出最佳 lambda
optimal_lambda <- lambda_seq[which.min(cv_errors)]

# 使用最佳 lambda 擬合模型
lasso_model <- coordinate_descent(X, y, optimal_lambda)

# 印出結果
print(lasso_model)