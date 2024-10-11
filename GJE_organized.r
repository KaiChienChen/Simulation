# function
row_switch  <- function(matrix1, r) {
  row1 <- matrix1[r, ]
  switched <- FALSE
  for (k in r:nrow(matrix1)) {
    if (matrix1[k, r] != 0) {
      matrix1[r, ] <- matrix1[k, ]
      matrix1[k, ] <- row1
      switched <- TRUE
      break
    }
  }
  if (!switched) {
    warning("No non-zero pivot found; matrix might be singular.")
  }
  return(matrix1)
}

gjrp <- function(matrix1,r){
  matrix1[r, ] <- matrix1[r, ] / matrix1[r, r]
  # r不是最後一行就進行row_process
  if (r!=nrow(matrix1)) {
    for (k in 1:(nrow(matrix1)-r)) {
      # 2. 把下面(k,i)都變成0
      matrix1[r+k,] <- matrix1[r+k,]- matrix1[r,]*matrix1[r+k,r]
      
    }
  }
  return(matrix1)
}

Gauss_Jordan <- function(matrix1) {
  for (i in 1:nrow(matrix1)) {
    if(matrix1[i,i]==0) {matrix1<- row_switch(matrix1,r=i)}
    matrix1 <- gjrp(matrix1 = matrix1,r = i)
  }
  return(matrix1)
}
Gauss_Jordan2 <- function(matrix1) {
  for (row1 in 2:nrow(matrix1)) {
    for (i in (row1-1):1) {
      matrix1[i,] <- matrix1[i,]-matrix1[row1,]*matrix1[i,row1] #往上減
    }
  }
  return(matrix1)
}

CreateIM <- function(matrix1){
  if(nrow(matrix1)==ncol(matrix1)){
      IdenM1 <- matrix(rep(0,nrow(matrix1)*ncol(matrix1)), nrow = nrow(matrix1))
  for(i in 1:nrow(IdenM1)){
    IdenM1[i,i] <- 1
  }

  }else {
     warning("Cannot Create a Identity Matrix")
  }
  return(IdenM1)
}

FindInver <- function(matrix1){
  Result1 <-Gauss_Jordan2(Gauss_Jordan(cbind(matrix1,CreateIM(matrix1))));Result1
  Inver <- (nrow(CreateIM(matrix1))+1):(2*nrow(CreateIM(matrix1)));Inver
  InverM <- Result1[,Inver]

  return(InverM)
}

# Test

m2 <- matrix(c(2,1,1,1,2,1,1,1,2),nrow = 3,ncol = 3,byrow = T);m2
FindInver(m2)
Gauss_Jordan(m2)
Gauss_Jordan2(Gauss_Jordan(m2))
m3 <- matrix(c(1,-2,3,9,-1,3,0,-4,2,-5,5,17), nrow = 3, ncol = 4,byrow = T);m3
Gauss_Jordan2(Gauss_Jordan(m3))

m4 <- matrix(c(1,3,1,9,1,1,-1,1,3,11,5,35), nrow = 3, ncol = 4,byrow = T);m4
Gauss_Jordan(m4)
FindInver(m4)
m4 <- gjrp(m4,2);m4
m5 <- gjrp(m4,2);m5
m5 <- row_switch(m4,3);m5
m5 <- gjrp(m5,2);m5

c1 <- cbind(m2,IdenM);c1
Gauss_Jordan(c1)
Gauss_Jordan2(m3)
m6 <- matrix(c(1,2,-1,6,4,2,4,-1,5),nrow = 3);m6
Gauss_Jordan(m6)
CreateIM(m6)
FindInver(m6) %*% m6
