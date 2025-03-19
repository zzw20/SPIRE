# We provide code for CC estimator under the normal setting with approximately 80% covariate censoring.
library(MASS)
library(nleqslv)
library(numDeriv)

for (i in 1:1000) {
  
  # date generation
  n <- 1000
  set.seed(i)
  z <- rbinom(n,1,0.5)
  c <- runif(n,min=z-0.5,max=z+0.5)
  
  ## generate X: X|C,Z follows a normal distribution
  diffx <- -1/2
  x <- rnorm(n, mean=c-diffx, sd=sqrt(z+1)/2)
  
  ## calculate the censoring rate
  cens.rate <- sum(c<=x)/n
  
  ## generate W and Delta
  w <- pmin(c,x)
  delta <- as.numeric(x<=c)
  
  ## Generate Y: Y=beta_0 + beta_1 x + beta_2 z + epsilon
  beta0 <- 0.5
  beta1 <- 0.2
  beta2 <- -0.2
  y <- beta0 + beta1*x + beta2*z + rnorm(n)
  beta <- c(beta0,beta1,beta2)
  
  # calculate CC estimator
  S.beta.f <- function(beta,y,x,z){
    
    res1 <- y-beta[1]-beta[2]*x-beta[3]*z
    res2 <- res1*x
    res3 <- res1*z
    c(res1,res2,res3)
    
  }
  cce <- function(beta,y,w,delta,z){
    res <- mapply(S.beta.f, y=y,x=w,z=z, MoreArgs = list(beta=beta))
    apply(res %*% diag(delta), MARGIN = 1,sum)
  }
  
  ## output CC estimator
  beta.hat <- nleqslv(beta,cce,y=y,w=w,delta=delta,z=z)$x
  
  # calculate estimated variance of CC estimator
  cce.func <- function(beta,y,w,delta,z){
    delta*S.beta.f(beta,y,w,z)
  }
  cce.deriv <- function(beta,y,w,delta,z){
    numDeriv::jacobian(function(beta) cce.func(beta,y,w,delta,z),beta)
  }
  res.cce <- mapply(cce.func, y=y,w=w,z=z,delta=delta,MoreArgs = list(beta=beta.hat))
  vn.hat <- res.cce %*% t(res.cce)/n
  jn.hat <-0
  for (i in 1:1000) {
    jn.hat <- jn.hat + cce.deriv(beta.hat,y[i],w[i],delta[i],z[i])
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))
  
  ## output estimated variance
  emp.cov/n
  
}
