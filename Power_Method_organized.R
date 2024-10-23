power_method <- function(matrix1, max_iter = 1000, tol = 1e-10) {
  # 任意選取向量 u1，並進行正規化
  u1 <- runif(nrow(matrix1))
  u1 <- u1 / sqrt(sum(u1^2))
  iter <- 0  # 初始化迭代計數
  eigenvalue0 <- 0  # 初始化特徵值
  differ <- tol + 1  # 初始化差異，使其大於容差
  while (iter <= max_iter && differ > tol) {
    iter <- iter + 1
    # 計算 A * u1
    result1 <- matrix1 %*% u1
    # 計算特徵值（向量的範數）
    eigenvalue1 <- sqrt(sum(result1^2))
    # 正規化向量 u1
    u1 <- result1 / eigenvalue1
    # 計算特徵值的變化量
    differ <- abs(eigenvalue1 - eigenvalue0)
    # 更新特徵值
    eigenvalue0 <- eigenvalue1
  }
  # 返回結果：最大特徵值、特徵向量和迭代次數
  return(list(eigenvalue = eigenvalue1, eigenvector = u1, iterations = iter))
}

# 測試矩陣
a <- matrix(c(2, 1, 1, 3), nrow = 2, byrow = TRUE)

# 執行 Power Method
result <- power_method(m1)

# 打印結果
cat("最大特徵值：", result$eigenvalue, "\n")
cat("對應的特徵向量：", result$eigenvector, "\n")
cat("迭代次數：", result$iterations, "\n")
