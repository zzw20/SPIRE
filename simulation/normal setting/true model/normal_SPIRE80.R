# We provide code for SPIRE under the normal setting 
# with approximately 80% covariate censoring and true working model.
library(MASS)
library(nleqslv)
library(numDeriv)
library(gaussquad)

for (i in 1:1000) {
  
  # date generation
  n <- 1000
  set.seed(i)
  z <- rbinom(n,1,0.5)
  c <- runif(n,min=z-0.5,max=z+0.5)
  
  ## generate X: X|C,Z follows a normal distribution
  diffx <- -0.5
  x <- rnorm(n, mean=c-diffx, sd=sqrt(z+1)/2)
  
  ## calculate the censoring rate
  cens.rate <- sum(c<=x)/n
  
  ## generate W and Delta
  w <- pmin(c,x)
  delta <- as.numeric(x<=c)
  
  ## set grid points for solving integral equation
  m <- 35
  x.a <- seq(from=mean(x)-3*sd(x), to=mean(x)+3*sd(x), length.out=m)
  
  ## Generate Y: Y=beta_0 + beta_1 x + beta_2 z + epsilon
  beta0 <- 0.5
  beta1 <- 0.2
  beta2 <- -0.2
  y <- beta0 + beta1*x + beta2*z + rnorm(n)
  beta <- c(beta0,beta1,beta2)
  
  # calculate SPIRE
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
  f.x.cz <- function(x,c,z){
    dnorm(x,mean=c-diffx, sd=1/2*sqrt(z+1))
  }
  S.beta <- function(beta,y,w,delta,z){
    if(delta==0){
      res1 <-integrate(function(x) S.beta.f1(beta,y,x,z)*f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),lower = w, upper = Inf
      )$value/integrate(function(x) f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),
                        lower=w, upper=Inf)$value
      res2 <-  integrate(function(x) S.beta.f2(beta,y,x,z)*f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),lower = w, upper = Inf
      )$value/integrate(function(x) f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),
                        lower=w, upper=Inf)$value
      res3 <- integrate(function(x) S.beta.f3(beta,y,x,z)*f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),lower = w, upper = Inf
      )$value/integrate(function(x) f.y.xz(beta,y,x,z)*f.x.cz(x,w,z),
                        lower=w, upper=Inf)$value
      c(res1,res2,res3)
    }else{
      S.beta.f(beta,y,w,z)
    }}
  b.cz <- function(beta,c,x.a,z){
    p <- dnorm(x.a, mean=c-diffx,sd=1/2*sqrt(z+1))/sum(dnorm(x.a, mean=c-diffx,sd=1/2*sqrt(z+1)))
    temp2.1 <- function(beta,y,c,x.a,z){
      sum((c<x.a)*S.beta.f1(beta,y,x.a,z)*f.y.xz(beta,y,x.a,z)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z)*p)
    }
    temp2.2 <- function(beta,y,c,x.a,z){
      sum((c<x.a)*S.beta.f2(beta,y,x.a,z)*f.y.xz(beta,y,x.a,z)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z)*p)
    }
    temp2.3 <- function(beta,y,c,x.a,z){
      sum((c<x.a)*S.beta.f3(beta,y,x.a,z)*f.y.xz(beta,y,x.a,z)*p)/sum((c<x.a)*f.y.xz(beta,y,x.a,z)*p)
    }
    cc <- gaussHermite(10)
    temp3 <- function(beta,c,x.a,x.aj,z){
      res1 <- sum(cc$w * sapply(cc$x, function(y) temp2.1(beta,sqrt(2)*y+beta[1]+beta[2]*x.aj+beta[3]*z,c,x.a,z)/sqrt(pi)))
      res2 <-  sum(cc$w * sapply(cc$x, function(y) temp2.2(beta,sqrt(2)*y+beta[1]+beta[2]*x.aj+beta[3]*z,c,x.a,z)/sqrt(pi)))
      res3 <- sum(cc$w * sapply(cc$x, function(y) temp2.3(beta,sqrt(2)*y+beta[1]+beta[2]*x.aj+beta[3]*z,c,x.a,z)/sqrt(pi)))
      c(res1,res2,res3)
    }
    mapply(x.aj=x.a,FUN=temp3,MoreArgs=list(c=c,beta=beta,z=z,x.a=x.a))
  }
  
  A.cz <- function(beta,c,x.a,z){
    p<-dnorm(x.a, mean=c-diffx,sd=1/2*sqrt(z+1))/sum(dnorm(x.a, mean=c-diffx,sd=1/2*sqrt(z+1)))
    temp4 <- function(y){
      (c<x.a)*f.y.xz(beta,y,x.a,z)*p/sum((c<x.a)*f.y.xz(beta,y,x.a,z)*p)
    }
    cc <- gaussHermite(10)
    mat.res <-  sapply(x.a,function(x.aj) {t(cc$w) %*% t(sapply(cc$x, function(y) temp4(sqrt(2)*y+beta[1]+beta[2]*x.aj+beta[3]*z)/sqrt(pi)))})
    t(mat.res)
  }
  a.0 <- function(beta,w,delta,x.a,z){
    if (delta==0){
      t(ginv(A.cz(beta,w,x.a,z)) %*% t(b.cz(beta,w,x.a,z)))
    }else{
      matrix(0,length(beta), m)
    }
  }
  
  spire.func <- function(beta,y,w,delta,z,x.a){
    if(delta==0){
      p<-dnorm(x.a, mean=w-diffx,sd=1/2*sqrt(z+1))/sum(dnorm(x.a, mean=c-diffx,sd=1/2*sqrt(z+1)))
      as.vector(S.beta(beta,y,w,delta,z)-a.0(beta,w,delta,x.a,z)%*%
                  ((w<x.a)*p*f.y.xz(beta,y,x.a,z))/sum((w<x.a)*p*f.y.xz(beta,y,x.a,z)))
    }else{
      S.beta(beta,y,w,delta,z)
    }
  }
  spire <- function(beta,y,w,delta,z,x.a){
    res <- mapply(spire.func,y=y,w=w,z=z,delta=delta,MoreArgs = list(beta=beta,x.a=x.a))
    apply(res,MARGIN = 1,sum)
  }
 
  ## output SPIRE
  beta.hat <- nleqslv::nleqslv(beta,spire,y=y,w=w,delta=delta,z=z,x.a=x.a)$x
  
  # calculate estimated variance of SPIRE
  spire.deriv <- function(beta,y,w,delta,z,x.a){
    numDeriv::jacobian(function(beta) spire.func(beta,y,w,delta,z,x.a),beta)
  }
  res.spire <- mapply(spire.func, y=y,w=w,z=z,delta=delta,MoreArgs = list(beta=beta.hat,x.a=x.a))
  vn.hat <- res.spire %*% t(res.spire)/n
  jn.hat <-0
  for (i in 1:1000) {
    jn.hat <- jn.hat + spire.deriv(beta.hat,y[i],w[i],delta[i],z[i],x.a)
    
  }
  jn.hat <- jn.hat/n
  emp.cov <- solve(jn.hat) %*% vn.hat %*%solve(t(jn.hat))
  
  ## output estimated variance
  emp.cov/n
}
