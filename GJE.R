# Overall test
## example 1
m1 <- matrix(data = c(0,1,1,1,2,1,1,1,2),nrow = 3,ncol = 3,byrow = T)

# row switch
## test
m1 <- matrix(data = c(2,1,1,1,0,1,1,1,2),nrow = 3,ncol = 3,byrow = T);m1
#################如果出現有(r,r)=0的情形，往下找非0的row(k,r)去做交換############
"row1 <- m1[1,]
for (k in 1:nrow(m1)) {
  if (m1[k,1]!=0) {m1[1,]<-m1[k,];
  m1[k,]<-row1;break}
    
}"
######### function ###########

row_switch  <- function(matrix1,r) {
  row1 <- matrix1[r,]
  for (k in r:nrow(matrix1)) {
    if (matrix1[k,r]!=0) {matrix1[r,]<-matrix1[k,];
    matrix1[k,]<-row1;break}
  }
  return(matrix1)
}
#############################

m1 <-row_switch(m1,r = 1)  ;m1 ##將r這個row做switch之動作if 符合以上說明之條件
  

# row process
## test
m1 <- matrix(data = c(2,1,1,1,0,1,1,1,2),nrow = 3,ncol = 3,byrow = T);m1
m1 <- gjrp(m1,r = 1);m1
m1 <- gjrp(m1,r = 2);m1
m1 <- gjrp(m1,r = 3);m1

# 1. 如果主對角元素(i,i)不為1，則將整行除以該元素使主對角元素變為1
m1[1, ] <- m1[1, ] / m1[1, 1];m1
# 2. 把下面(k,i)都變成0
m1[1+2,] <- m1[1+2,]- m1[1,]*m1[1+2,1];m1
## for loop
r1<- 1
m1[r1, ] <- m1[r1, ] / m1[r1, r1];m1
for (k in 1:(nrow(m1)-r1)) {
  # 2. 把下面(k,i)都變成0
  m1[r1+k,] <- m1[r1+k,]- m1[r1,]*m1[r1+k,r1];m1
  
}
m1
########### function ##############

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
#############################

# Gauss-Jordan Elimination

###### function ###############

Gauss_Jordan <- function(matrix1) {
  for (i in 1:nrow(matrix1)) {
    if(matrix1[i,i]==0) {matrix1<- row_switch(matrix1,r=i)}
    matrix1 <- gjrp(matrix1 = matrix1,r = i)
  }
  return(matrix1)
}
###############################
for (i in 1:nrow(m1)) {
  if(m1[i,i]==0) {m1<- row_switch(m1,r=i)}
  m1 <- gjrp(matrix1 = m1,r = i)
}
m1

# Gauss-Jordan Elimination 2 轉換成identity matrix
## example
c1 <- Gauss_Jordan(m2);c1
c1[3,3]
for (row1 in 2:nrow(c1)) {
  for (i in (row1-1):1) {
    c1[i,] <- c1[i,]-c1[row1,]*c1[i,row1];c1#往上減
  }
}
c1
########### function ############
Gauss_Jordan2 <- function(matrix1) {
  for (row1 in 2:nrow(matrix1)) {
    for (i in (row1-1):1) {
      matrix1[i,] <- matrix1[i,]-matrix1[row1,]*matrix1[i,row1] #往上減
    }
  }
  return(matrix1)
}
#################################


########## test ################
IdenM <- matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3);IdenM

m1 <- matrix(c(2,1,-1,8,-3,-1,2,-11,-2,1,2,-3),nrow = 3,ncol = 4,byrow = T);m1
m2 <- matrix(c(2,1,1,1,2,1,1,1,2),nrow = 3,ncol = 3,byrow = T);m2

m3 <- matrix(c(1,-2,3,9,-1,3,0,-4,2,-5,5,17), nrow = 3, ncol = 4,byrow = T);m3
Gauss_Jordan2(Gauss_Jordan(m3))

m4 <- matrix(c(1,3,1,9,1,1,-1,1,3,11,5,35), nrow = 3, ncol = 4,byrow = T);m4
Gauss_Jordan(m4)
m4 <- gjrp(m4,2);m4
m5 <- gjrp(m4,2);m5
m5 <- row_switch(m4,3);m5
m5 <- gjrp(m5,2);m5

c1 <- cbind(m2,IdenM);c1
Gauss_Jordan(c1)
Gauss_Jordan2(m3)
################### adjusted row_switch function ###############
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
