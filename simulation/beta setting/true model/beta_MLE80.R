# We provide code for MLE under the beta setting 
# with approximately 85% covariate censoring and true working model.
library(MASS)
library(nleqslv)
library(numDeriv)

for (i in 1:1000) {
  
  # date generation
  n <- 3000
  set.seed(i)
  z1.par <- c(3.5383,11.4963)
  z1 <- rbeta(n,z1.par[1],z1.par[2])
  z2 <- rbinom(n,1,0.5)
  alpha.c <- 0.3+z1
  beta.c <- 1.1+z2
  c <- rbeta(n,alpha.c,beta.c)
  
  ## generate X: X|C,Z follows a beta distribution
  alpha.x <- 1.6+5*c
  beta.x <- 2+z1+z2
  x <- rbeta(n,alpha.x,beta.x)
  
  ## calculate the censoring rate
  cens.rate <- sum(c<=x)/n
  
  ## generate W and Delta
  w <- pmin(c,x)
  delta <- as.numeric(x<=c)
  
  ## Generate Y: Y=beta_0 + beta_1 x + beta_2 z1 + beta_3 z2 + beta_4 xz_2 + epsilon
  beta0 <- 1.3
  beta1 <- -1.8
  beta2 <- -1.5
  beta3 <- 0.1
  beta4 <- 0.2
  y <- beta0 + beta1*x + beta2*z1 + beta3*z2+beta4*x*z2+rnorm(n,sd=1)
  beta <- c(beta0,beta1,beta2,beta3,beta4,1)
  
  # calculate MLE
  S.beta.f <- function(beta,y,x,z1,z2){
    res1 <- (y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
    res2 <- res1*x
    res3 <- res1*z1
    res4 <- res1*z2
    res5 <- res1*x*z2
    res6 <- -1/(2*beta[6])+(y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)^2/(2*beta[6]^2)
    c(res1,res2,res3,res4,res5,res6)
  }
  S.beta.f1 <- function(beta,y,x,z1,z2){
    (y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
  }
  S.beta.f2 <- function(beta,y,x,z1,z2){
    eta <-(y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
    eta*x
  }
  S.beta.f3 <- function(beta,y,x,z1,z2){
    eta <-(y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
    eta*z1
  }
  S.beta.f4 <- function(beta,y,x,z1,z2){
    eta <- (y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
    eta*z2
  }
  S.beta.f5 <- function(beta,y,x,z1,z2){
    eta <-(y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)/beta[6]
    eta*x*z2
  }
  S.beta.f6 <- function(beta,y,x,z1,z2){
    -1/(2*beta[6])+(y-beta[1]-beta[2]*x-beta[3]*z1-beta[4]*z2-beta[5]*x*z2)^2/(2*beta[6]^2)
  }
  f.x.cz <- function(x,c,z1,z2){
    dbeta(x,1.6+5*c,2+z1+z2)
  }
  f.y.xz <- function(beta,y,x,z1,z2){
    eta <- beta[1]+beta[2]*x+beta[3]*z1+beta[4]*z2+beta[5]*x*z2
    dnorm(y, mean = eta, sd=sqrt(beta[6]))
  }
  mle.func <- function(beta,y,w,delta,z1,z2){
    if(delta==0){
      res1 <-integrate(function(x) S.beta.f1(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      res2 <-integrate(function(x) S.beta.f2(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      res3 <- integrate(function(x) S.beta.f3(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      res4 <- integrate(function(x) S.beta.f4(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      res5 <- integrate(function(x) S.beta.f5(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      res6 <- integrate(function(x) S.beta.f6(beta,y,x,z1,z2)*f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),lower = w, upper = 1
      )$value/integrate(function(x) f.y.xz(beta,y,x,z1,z2)*f.x.cz(x,w,z1,z2),
                        lower=w, upper=1)$value
      c(res1,res2,res3,res4,res5,res6)
    }else{
      S.beta.f(beta,y,w,z1,z2)
    }}
  
  mle<- function(beta,y,w,delta,z1,z2){
    res <- mapply(mle.func, y=y,w=w,z1=z1,z2=z2, delta=delta,MoreArgs = list(beta=beta))
    apply(res,MARGIN = 1,sum)
  }
  
  ## output MLE
  beta.hat <- nleqslv::nleqslv(beta,mle,y=y,w=w,delta=delta,z1=z1,z2=z2)$x
  
  # calculate estimated variance of MLE
  mle.deriv <- function(beta,y,w,delta,z1,z2){
    numDeriv::jacobian(function(beta) mle.func(beta,y,w,delta,z1,z2),beta)
  }
  res.mle <- mapply(mle.func, y=y,w=w,z1=z1,z2=z2,delta=delta,MoreArgs = list(beta=beta.hat))
  vn.hat <- res.mle %*% t(res.mle)/n
  jn.hat <-0
  for (i in 1:3000) {
    jn.hat <- jn.hat + mle.deriv(beta.hat,y[i],w[i],delta[i],z1[i],z2[i])
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))
  
  ## output estimated variance
  emp.cov/n
}
