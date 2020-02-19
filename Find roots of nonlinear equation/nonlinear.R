FixedPoint <- function(f, start, tol = 1e-8, max.iter = 20){
  iter <- 0; error <- tol * 2 
  while(error > tol){
    temp <- start
    start <- f(temp)
    error <- abs(start - temp)
    iter <- iter + 1
    if(iter > max.iter){
      stop("Does not converge!")
    }
  }
  return(start)
}

Bisection <- function(f, x1, x2, tol = 1e-8){
  f1 <- f(x1); f2 <- f(x2)
  if(abs(f1) < tol) return(x1)
  if(abs(f2) < tol) return(x2)
  if(sign(f1) * sign(f2) > 0){
    stop("The initial condition is not satisfied")
  }
  error <- tol * 2
  while(error > tol){
    x3 <- (x1 + x2) / 2
    f3 <- f(x3)
    if(abs(f3) < tol) return(x3)
    if(sign(f2) * sign(f3) > 0){
      error <- x2 - x1
      x2 <- x3; f2 <- f3
    }else{
      error <- x2 - x1
      x1 <- x3; f1 <- f3
    }
    print(error)
  }
  root <- (x1 + x2) / 2
  return(root)
}

Newton <- function(f, df, start, tol = 1e-8, max.iter = 20, bound = 50){
  iter <- 0; error <- tol * 2
  while(error > tol){
    temp <- start
    dif <- f(start) / df(start)
    start <- start - dif
    error <- abs(dif) / abs(temp)
    iter <- iter + 1
    if((abs(start) > bound) | (iter > max.iter)){
      stop("Does not converge!")
    }
  }
  return(start)
}

Newton2 <- function(f, df, ddf, start, tol = 1e-8, max.iter = 20, bound = 50){
  iter <- 0; error <- tol * 2
  while(error > tol){
    temp <- start
    dif <- f(start) * df(start) / (df(start)^2 - f(start) * ddf(start))
    start <- start - dif
    error <- abs(dif)
    iter <- iter + 1
    if((abs(start) > bound) | (iter > max.iter)){
      stop("Does not converge!")
    }
  }
  return(start) 
}

NewtonRhapson <- function(f, df, start, tol = 1e-8, ddf = NA, second = FALSE,
                          max.iter = 20, bound = 50){
  if(!second){
    Newton(f, df, start, tol, max.iter, bound)
  }else{
    if(is.na(ddf)){
      stop("Need sedond derivative")
    }else{
      Newton2(f, df, ddf, start, tol, max.iter, bound)
    }
  }
}  
  
Secant <- function(f, x1, x2, tol = 1e-8, max.iter = 40, bound = 50){
  iter <- 0; error <- tol * 2
  while(error > tol){
    secant <- (f(x2) - f(x1)) / (x2 - x1)
    dif <- f(x2) / secant
    x1 < - x2; x2 <- x2 - dif
    error <- abs(dif)
    iter <- iter + 1
    if((abs(x2) > bound) | (iter > max.iter)){
      stop("Does not converge!")
    }
  }
  return(st2)
}

Brent <- function(f, a, b, tol = 1e-8, max.iter = 30){
  x1 <- a; x2 <- b; iter <- 0
  f1 <- f(x1);  f2 <- f(x2)
  if(abs(f1) < tol) return(x1)
  if(abs(f2) < tol) return(x2)
  if(sign(f1) * sign(f2) > 0){
    stop("The initial condition is not satisfied")
  }
  x3 <- (x1 + x2) / 2
  while(iter < max.iter){
    f3 <- f(x3)
    if(abs(f3) < tol) return(x3)
    if(sign(f1) * sign(f3) < 0){
      b <- x3
    }else{
      a <- x3
    }
    if((b - a) < tol * max(abs(b), 1)) return((a + b)/2)
    denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
    numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3) + f1 * x2 * (f3 - f1)
    logi = (denom == 0)
    dx <- f3 * numer / denom *(1- logi) + (b - a) * logi 
    x <- x3 + dx
    if((b - x) * (x - a) < 0){
      dx <- (b - a)/2
      x <- a + dx
    }
    if(x < x3){
      x2 <- x3; f2 <- f3
    }else{
      x2 <- x3; f1 <- f3
    }
    x3 <- x
    iter <- iter + 1
  }
  stop("Does not converge!")
}