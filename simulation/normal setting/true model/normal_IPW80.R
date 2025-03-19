# We provide code for IPW estimator under the normal setting
# with approximately 80% covariate censoring and true working model.
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
  
  # calculate IPW estimator
  S.beta.f <- function(beta,y,x,z){
    
    res1 <- y-beta[1]-beta[2]*x-beta[3]*z
    res2 <- res1*x
    res3 <- res1*z
    c(res1,res2,res3)
    
  }
  ipe.func <- function(beta,y,w,delta,z){
    if(delta==0){
      c(0,0,0)
    }else{
      lo <- max(z-0.5,w)
      norm.xz <- pnorm(z+0.5, mean=w+diffx,sd=sqrt(z+1)/2)-pnorm(z-0.5, mean=w+diffx,sd=sqrt(z+1)/2)
      pr <- (pnorm(z+0.5, mean=w+diffx,sd=sqrt(z+1)/2)-pnorm(lo, mean=w+diffx,sd=sqrt(z+1)/2))/norm.xz
      S.beta.f(beta,y,w,z)/pr
    }
  }
  ipe <- function(beta,y,w,delta,z){
    res <- mapply(ipe.func, y=y,w=w,z=z,delta=delta, MoreArgs = list(beta=beta))
    apply(res,MARGIN = 1,sum)
  }
  
  ## output IPW estimator
  beta.hat <- nleqslv(beta,ipe,y=y,w=w,delta=delta,z=z)$x
  
  # calculate estimated variance of IPW estimator
  ipe.deriv <- function(beta,y,w,delta,z){
    numDeriv::jacobian(function(beta) ipe.func(beta,y,w,delta,z),beta)
  }
  res.ipe <- mapply(ipe.func, y=y,w=w,z=z,delta=delta,MoreArgs = list(beta=beta.hat))
  vn.hat <- res.ipe %*% t(res.ipe)/n
  jn.hat <-0
  for (i in 1:1000) {
    jn.hat <- jn.hat + ipe.deriv(beta.hat,y[i],w[i],delta[i],z[i])
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))

  ## output estimated variance
  emp.cov/n
}
