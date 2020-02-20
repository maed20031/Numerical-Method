dot <- function(a, b){
  return(as.numeric(t(a) %*% b))
}

GaussSeidel <- function(A, x, b, tol = 1e-8){
  L <- U <- A
  L[upper.tri(A, TRUE)]<- 0
  U[lower.tri(A, TRUE)]<- 0
  D <- diag(diag(A)^(-1))
  while(TRUE){
    x_new <- D %*%(b - (L+U)%*% x)
    err <- x_new - x
    if(dot(err, err) < tol){
      break
    }
    x <- x_new
  }
  return(x_new)
}

ConGrad <- function(A, x, b, tol = 1e-8){
 n <- length(b); r <- b - A %*% x; s <- r
 for(i in 1:n){
   u <- A %*% s
   alpha <- dot(s, r) / dot(s, u)
   x <- x + alpha * s
   r <- b - A %*% x
   if(sqrt(dot(r, r)) < tol){
     break
   }else{
     beta <- -dot(r, u) / dot(s, u)
     s <- r + beta * s
   }
 }
 return(x)
}

A <- matrix(c(4, -1, 1, -1, 4, -2, 1, -2, 4),3)
x <- c(0,0,0)
b <- c(12, -1, 5)
