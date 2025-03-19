# We provide code for SPIRE under the beta setting 
# with approximately 85% covariate censoring and uniform working model.
library(MASS)
library(nleqslv)
library(numDeriv)
library(gaussquad)

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
  
  ## set grid points for solving integral equation
  m <- 20
  x.a <- seq(from=min(w), to=0.9999, length.out=m)
  p <- rep(1/m,m)
  
  ## Generate Y: Y=beta_0 + beta_1 x + beta_2 z1 + beta_3 z2 + beta_4 xz_2 + epsilon
  beta0 <- 1.3
  beta1 <- -1.8
  beta2 <- -1.5
  beta3 <- 0.1
  beta4 <- 0.2
  y <- beta0 + beta1*x + beta2*z1 + beta3*z2+beta4*x*z2+rnorm(n,sd=1)
  beta <- c(beta0,beta1,beta2,beta3,beta4,1)
  
  # calculate SPIRE
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
  S.beta <- function(beta,y,w,delta,z1,z2){
    if(delta==0){
      denom <- sum(f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)
      res1 <- sum(S.beta.f1(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      res2 <- sum(S.beta.f2(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      res3 <- sum(S.beta.f3(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      res4 <- sum(S.beta.f4(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      res5 <- sum(S.beta.f5(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      res6 <- sum(S.beta.f6(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*(x.a>w)*p)/denom
      c(res1,res2,res3,res4,res5,res6)
    }else{
      S.beta.f(beta,y,w,z1,z2)
    }}
  b.cz <- function(beta,c,x.a,z1,z2){
    temp2.1 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f1(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    temp2.2 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f2(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    temp2.3 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f3(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    temp2.4 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f4(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    temp2.5 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f5(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    temp2.6 <- function(beta,y,c,x.a,z1,z2){
      sum((c<x.a)*S.beta.f6(beta,y,x.a,z1,z2)*f.y.xz(beta,y,x.a,z1,z2)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    cc <- gaussHermite(10)
    temp3 <- function(beta,c,x.a,x.aj,z1,z2){
      res1 <- sum(cc$w * sapply(cc$x, function(y) temp2.1(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      res2 <-  sum(cc$w * sapply(cc$x, function(y) temp2.2(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      res3 <- sum(cc$w * sapply(cc$x, function(y) temp2.3(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      res4 <- sum(cc$w * sapply(cc$x, function(y) temp2.4(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      res5 <- sum(cc$w * sapply(cc$x, function(y) temp2.5(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      res6 <- sum(cc$w * sapply(cc$x, function(y) temp2.6(beta,sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2,c,x.a,z1,z2)/sqrt(pi)))
      c(res1,res2,res3,res4,res5,res6)
    }
    mapply(x.aj=x.a,FUN=temp3,MoreArgs=list(c=c,beta=beta,z1=z1,z2=z2,x.a=x.a))
  }
  A.cz <- function(beta,c,x.a,z1,z2){
    temp4 <- function(y){
      (c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    cc <- gaussHermite(10)
    mat.res <-  sapply(x.a,function(x.aj) {t(cc$w) %*% t(sapply(cc$x, function(y) temp4(sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2)/sqrt(pi)))})
    t(mat.res)
  }
  a.0 <- function(beta,w,delta,x.a,z1,z2){
    if (delta==0){
      t(ginv(A.cz(beta,w,x.a,z1,z2)) %*% t(b.cz(beta,w,x.a,z1,z2)))
    }else{
      matrix(0,length(beta), m)
    }
  }
  spire.func <- function(beta,y,w,delta,z1,z2,x.a){
    if(delta==0){
      as.vector(S.beta(beta,y,w,delta,z1,z2)-a.0(beta,w,delta,x.a,z1,z2)%*%
                  ((w<x.a)*p*f.y.xz(beta,y,x.a,z1,z2))/sum((w<x.a)*p*f.y.xz(beta,y,x.a,z1,z2)))
    }else{
      S.beta(beta,y,w,delta,z1,z2)
    }
  }
  spire <- function(beta,y,w,delta,z1,z2,x.a){
    res <- mapply(spire.func, y=y,w=w,z1=z1,z2=z2,delta=delta,MoreArgs = list(beta=beta,x.a=x.a))
    apply(res,MARGIN = 1,sum)
  }
  
  ## output SPIRE
  beta.hat <- nleqslv::nleqslv(beta,spire,y=y,w=w,delta=delta,z1=z1,z2=z2,x.a=x.a)$x
  
  # calculate estimated variance of SPIRE
  spire.deriv <- function(beta,y,w,delta,z1,z2,x.a){
    numDeriv::jacobian(function(beta) spire.func(beta,y,w,delta,z1,z2,x.a),beta)
  }
  res.spire <- mapply(spire.func, y=y,w=w,z1=z1,z2=z2,delta=delta,MoreArgs = list(beta=beta.hat,x.a=x.a))
  vn.hat <- res.spire %*% t(res.spire)/n
  jn.hat <-0
  for (i in 1:3000) {
    jn.hat <- jn.hat + spire.deriv(beta.hat,y[i],w[i],delta[i],z1[i],z2[i],x.a)
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))
  
  ## output estimated variance
  emp.cov/n
}
