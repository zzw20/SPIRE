# We provide code for SPIRE under the beta setting 
# with approximately 85% covariate censoring and Kaplan-Meier working model.
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
  
  ## set grid points 
  index.o <- order(w)
  y.o <- y[index.o]
  w.o <- w[index.o]
  delta.o <- delta[index.o]
  w.o1 <- w.o[delta.o==1]
  z1.o <- z1[index.o]
  z2.o <- z2[index.o]
  m <- 20
  l <- 20
  x.a <- seq(from=0,to=1,length.out=m)
  z1.a <- seq(from=0,to=1,length.out=l)

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
  f.y.xz <- function(beta,y,x,z1,z2){
    eta <- beta[1]+beta[2]*x+beta[3]*z1+beta[4]*z2+beta[5]*x*z2
    dnorm(y, mean = eta, sd=sqrt(beta[6]))
  }
  kern = function(cond, cond_data, h){ 
    exp(- 0.5 * apply(as.matrix(cond_data), 1, function(x) sum((cond - x)^2)) / h^2) / (sqrt(2*pi)*h)
  }
  surv.x2 = function(t, w_data, delta_data,
                     z1, z2, z1_data, z2_data, h1){
    n = length(w_data)
    idx_j = which((w_data[z2_data == z2] <= t) & (delta_data[z2_data == z2] == 1)) #index of j in z_data==z
    if(length(idx_j)!=0){
      kernel_vals = kern(z1, z1_data[z2_data==z2], h1) 
      denom = sapply((w_data[z2_data == z2])[idx_j], function(w) {
        sum(kernel_vals * (w_data[z2_data==z2] >= w))})
      log_vals = log(pmax(1 - kernel_vals[idx_j] / denom, 1/n))
      exp(sum(log_vals))
    } else{
      1
    }
  }
  h<-0.05
  f.xz.hat <- function(z1.i,z2.i){
    temp.surv <- sapply(w.o1, function(t) surv.x2(t,w,delta,z1.i,z2.i,z1,z2,h))
    temp.dist <- 1-temp.surv
    temp.dist-c(0,temp.dist[1:length(temp.dist)-1])
  }
  z2.a <- unique(z2)
  f.xz.temp2 <- c()
  for (i in 1:length(z1.a)) {
    for (j in 1:length(z2.a)) {
      f.xz.temp2 <- cbind(f.xz.temp2,f.xz.hat(z1.a[i],z2.a[j]))
    }
  }
  f.xz.array <- array(f.xz.temp2, dim=c(length(w.o1),2,l)) # f.xz.hat(z1.u[i],z2.u[j])=f.xz.array[,j,i]
  S.beta <- function(beta,y,w,delta,z1,z2,p){
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
  b.cz <- function(beta,c,x.a,z1,z2,p){
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
  A.cz <- function(beta,c,x.a,z1,z2,p){
    temp4 <- function(y){
      (c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p/sum((c<x.a)*f.y.xz(beta,y,x.a,z1,z2)*p)
    }
    cc <- gaussHermite(10)
    mat.res <-  sapply(x.a,function(x.aj) {t(cc$w) %*% t(sapply(cc$x, function(y) temp4(sqrt(2*beta[6])*y+beta[1]+beta[2]*x.aj+beta[3]*z1+beta[4]*z2+beta[5]*x.aj*z2)/sqrt(pi)))})
    t(mat.res)
  }
  a.0. <- function(beta,w,delta,x.a,z1,z2,p){
    if (delta==0){
      t(ginv(A.cz(beta,w,x.a,z1,z2,p)) %*% t(b.cz(beta,w,x.a,z1,z2,p)))
    }else{
      matrix(0,length(beta), m)
    }
  }
  kern2 = function(cond, cond_data, h, weight){ 
    sum(exp(- 0.5 * apply(as.matrix(cond_data), 1, function(x) sum((cond - x)^2)) / h^2) / (sqrt(2*pi)*h)*weight)
  }
  
  h2 <- 0.05
  spire.func <- function(beta,y,w,delta,z1,z2,x.a){
    if(delta==0){
      i <- findintervals(z1,z1.a)
      j <- which(unique(z2.a)==z2)
      temp.p1 <- sapply(x.a, function(x) kern2(x,w.o1,h2,f.xz.array[,j,i]))
      temp.p2 <- sapply(x.a, function(x) kern2(x,w.o1,h2,f.xz.array[,j,i+1]))
      p <- temp.p1+(z1-z1.a[i])*(temp.p2-temp.p1)/(z1.a[i+1]-z1.a[i])
      if(sum((w<x.a)*p)>0){
        S.beta(beta,y,w,delta,z1,z2,p)-a.0.gauss3(beta,w,delta,x.a,z1,z2,p)%*%
          ((w <x.a)*p*f.y.xz(beta,y,x.a,z1,z2))/sum((w<x.a)*p*f.y.xz(beta,y,x.a,z1,z2))
      }else{
        c(0,0,0,0,0,0)
      }
    }
    else{
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
