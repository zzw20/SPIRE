# We provide code for MLE under the normal setting 
# with approximately 80% covariate censoring and uniform working model.
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
  
  # set grid points
  m <- 30
  p <- rep(1/m,m)
  x.a <- seq(from=mean(x)-3*sd(x), to=mean(x)+3*sd(x), length.out=m)
  
  ## Generate Y: Y=beta_0 + beta_1 x + beta_2 z + epsilon
  beta0 <- 0.5
  beta1 <- 0.2
  beta2 <- -0.2
  y <- beta0 + beta1*x + beta2*z + rnorm(n)
  beta <- c(beta0,beta1,beta2)
  
  # calculate MLE
  S.beta.f <- function(beta,y,x,z){
    res1 <- y-beta[1]-beta[2]*x-beta[3]*z
    res2 <- res1*x
    res3 <- res1*z
    c(res1,res2,res3)
  }
  S.beta.f1 <- function(beta,y,x,z){
    y-beta[1]-beta[2]*x-beta[3]*z
  }
  S.beta.f2 <- function(beta,y,x,z){
    eta <-y-beta[1]-beta[2]*x-beta[3]*z
    eta*x
  }
  S.beta.f3 <- function(beta,y,x,z){
    eta <-y-beta[1]-beta[2]*x-beta[3]*z
    eta*z
  }
  f.y.xz <- function(beta,y,x,z){
    eta <- beta[1]+beta[2]*x+beta[3]*z
    dnorm(y, mean = eta, sd=1)
  }
  
  mle.func <- function(beta,y,w,delta,z){
    if(delta==0){
      pr.c <- sum(mapply(function(x,p) f.y.xz(beta,y,x,z)*(x>w)*p,x=x.a,p=p))
      res1 <- sum(mapply(function(x,p) f.y.xz(beta,y,x,z)*S.beta.f1(beta,y,x,z)*(x>w)*p,x=x.a,p=p))/pr.c
      res2 <- sum(mapply(function(x,p) f.y.xz(beta,y,x,z)*S.beta.f2(beta,y,x,z)*(x>w)*p,x=x.a,p=p))/pr.c
      res3 <- sum(mapply(function(x,p) f.y.xz(beta,y,x,z)*S.beta.f3(beta,y,x,z)*(x>w)*p,x=x.a,p=p))/pr.c
      c(res1,res2,res3)
    }else{
      S.beta.f(beta,y,w,z)
    }
  } 
  mle<- function(beta,y,w,delta,z){
    res <- mapply(mle.func, y=y,w=w,z=z, delta=delta,MoreArgs = list(beta=beta))
    apply(res,MARGIN = 1,sum)
  }
  
  ## output MLE
  beta.hat <- nleqslv::nleqslv(beta,mle,y=y,w=w,delta=delta,z=z)$x
  
  # calculate estimated variance of MLE
  mle.deriv <- function(beta,y,w,delta,z){
    numDeriv::jacobian(function(beta) mle.func(beta,y,w,delta,z),beta)
  }
  res.mle <- mapply(mle.func, y=y,w=w,z=z,delta=delta,MoreArgs = list(beta=beta.hat))
  vn.hat <- res.mle %*% t(res.mle)/n
  jn.hat <-0
  for (i in 1:1000) {
    jn.hat <- jn.hat + mle.deriv(beta.hat,y[i],w[i],delta[i],z[i])
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))
  
  ## output estimated variance
  emp.cov/n
}
