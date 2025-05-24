rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select
hessian <- numDeriv::hessian

load(file="_R_Data_Files/full_dataset.RData")

# Split into training and testing
set.seed(9999) # set random seed for reproducibility
keep.half <- sample(1:nrow(dat),ceiling(nrow(dat)/2))
save.half <- !(1:nrow(dat) %in% keep.half)
dat.train <- dat[keep.half,]

# Perform model fits
## 1. One Theta; No Time
pdf.under.alt.symbolic <- function(x,theta){
  x1 <- x/2
  x2 <- 1-x/2
  1/2*exp(-1/2*theta*(theta + 2*sqrt(2)*erfinv(2*x1-1)))+
    1/2*exp(-1/2*theta*(theta + 2*sqrt(2)*erfinv(2*x2-1)))
}

fit.one.theta.no.time <- fitdistr(dat.train$p,
                pdf.under.alt.symbolic,
                start=list(theta=3.1),
                method="L-BFGS-B",
                lower=1,upper=5)

theta <- fit.one.theta.no.time$estimate[["theta"]]
ci <- confint(fit.one.theta.no.time)

one.theta.nt <- list(pars=theta,lower=ci[1],upper=ci[2],
                     loglik=fit.one.theta.no.time$loglik)

### 2. One Theta; Linear Time
pdf.diff.by.year <- function(x,time.diff,slope,theta0){
  theta <- theta0+slope*time.diff
  pdf.under.alt.symbolic(x,theta)
}

nll.fxn.by.year <- function(pars,dat.train){
  slope=pars[["slope"]]
  theta0=pars[["theta0"]]
  total.log.like <- 0
  for(ii in 1:nrow(dat.train)){
    tmp.like <- pdf.diff.by.year(
      x=dat.train[ii,"p"],
      time=dat.train[ii,"time.diff"],
      slope,theta0)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

par0 <- c(slope=0,
          theta0=fit.one.theta.no.time$estimate[["theta"]])
nll.fxn.by.year(pars=par0,dat.train=dat.train)

optim.out <- optim(par=par0,
                   nll.fxn.by.year,
                   dat.train=dat.train,
                   method="L-BFGS-B",
                   lower=c(-1,2),upper=c(1,4),
                   hessian=T)
pars <- optim.out$par
hess <- hessian(nll.fxn.by.year,optim.out$par,dat.train=dat.train)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- pars - z.crit*sigma
upper <- pars + z.crit*sigma

one.theta.lt <- list(pars=pars,lower=lower,upper=upper,
                     loglik=-optim.out$value)

## 3. Two-Dirac Delta Theta; No Time
two.hypth <- function(x,theta,null.weight){
  theta_star <- theta/(1-null.weight)
  null.weight*pdf.under.alt.symbolic(x,0)+
    (1-null.weight)*pdf.under.alt.symbolic(x,theta_star)
}
two.hypth.fit <- fitdistr(dat.train$p,
                          two.hypth,
                          start=list(theta=3,null.weight=0.15),
                          method="L-BFGS-B",
                          lower=c(1,0),upper=c(5,1))
pars <- two.hypth.fit$estimate
ci <- confint(two.hypth.fit)
two.theta.nt <- list(pars=pars,lower=ci[,1],upper=ci[,2],
                     loglik=two.hypth.fit$loglik)

## 4. Two-Dirac Delta Theta; Linear Time
pdf.diff.by.year.2d <- function(x,time.diff,slope,theta0,null.weight){
  theta <- theta0+slope*time.diff
  two.hypth(x,theta,null.weight)
}

nll.fxn.by.year <- function(pars,dat.train){
  slope=pars[["slope"]]
  theta0=pars[["theta0"]]
  null.weight=pars[["null.weight"]]
  total.log.like <- 0
  for(ii in 1:nrow(dat.train)){
    tmp.like <- pdf.diff.by.year.2d(
      x=dat.train[ii,"p"],
      time=dat.train[ii,"time.diff"],
      slope,theta0,null.weight)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

par0 <- c(slope=0,
          theta0=fit.one.theta.no.time$estimate[["theta"]],
          null.weight=two.hypth.fit$estimate[["null.weight"]])
nll.fxn.by.year(pars=par0,dat.train=dat.train)

optim.out <- optim(par=par0,
                   nll.fxn.by.year,
                   dat.train=dat.train,
                   method="L-BFGS-B",
                   lower=c(-1,2,0),upper=c(1,4,1))
pars <- optim.out$par
hess <- hessian(nll.fxn.by.year,optim.out$par,dat.train=dat.train)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- pars - z.crit*sigma
upper <- pars + z.crit*sigma

two.theta.lt <- list(pars=pars,lower=lower,upper=upper,
                     loglik=-optim.out$value)


## 5. Exponential Mixture Distribution 
int.fxn.exp <- function(theta,x,rate){
  pdf.under.alt.symbolic(x,theta)*dexp(theta,rate=rate)
}
like.fxn.one.exp <- function(x,rate){
  integrate(int.fxn.exp,0,Inf,x=x,rate=rate)$value
}
nll.fxn.exp <- function(pars,dat.train){
  theta=pars[["theta"]]
  rate <- 1/theta
  total.log.like <- 0
  for(ii in 1:length(dat.train)){
    tmp.like <- like.fxn.one.exp(dat.train[[ii]],rate=rate)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

par0 <- c(theta=2.8)
nll.fxn.exp(pars=par0,dat.train=dat.train$p)

optim.out <- optim(par=par0,
                   nll.fxn.exp,
                   dat.train=dat.train$p,
                   method="L-BFGS-B",
                   lower=c(2),upper=c(4))
est.rate <- optim.out$par
hess <- hessian(nll.fxn.exp,optim.out$par,dat.train=dat.train$p)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- est.rate - z.crit*sigma
upper <- est.rate + z.crit*sigma

exp.fit <- list(pars=est.rate,lower=lower,upper=upper)
exp.fit.rate <- lapply(exp.fit,function(x)1/x)
exp.fit$loglik <- -optim.out$value

## 6. Exponential Mixture Distribution, Linear Time 
int.fxn.exp.lt <- function(theta,x,time.diff,theta0,slope){
  rate <- 1/(theta0+slope*time.diff)
  pdf.under.alt.symbolic(x,theta)*dexp(theta,rate=rate)
}
like.fxn.one.exp.lt <- function(x,time.diff,theta0,slope){
  rate <- 1/(theta0+slope*time.diff)
  integrate(int.fxn.exp.lt,0,Inf,x=x,
            time.diff=time.diff,
            theta0=theta0,
            slope=slope)$value
}
nll.fxn.exp.lt <- function(pars,dat.train){
  theta0=pars[["theta0"]]
  slope=pars[["slope"]]
  total.log.like <- 0
  for(ii in 1:nrow(dat.train)){
    tmp.like <- like.fxn.one.exp.lt(dat.train$p[ii],dat.train$time.diff[ii],
                             theta0=theta0,slope=slope)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

t0 <- as.numeric(exp.fit$par)
par0 <- c(theta0=t0,slope=0.04)
nll.fxn.exp.lt(pars=par0,dat.train=dat.train)

optim.out <- optim(par=par0,
                   nll.fxn.exp.lt,
                   dat.train=dat.train,
                   method="L-BFGS-B",
                   lower=c(2,0),upper=c(4,0.05))
est.rate <- optim.out$par
hess <- hessian(nll.fxn.exp.lt,optim.out$par,dat.train=dat.train)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
sigma[2] <- sigma[2]/100
z.crit <- qnorm(0.975)
lower <- est.rate - z.crit*sigma
upper <- est.rate + z.crit*sigma

exp.fit.mean.lt <- list(pars=est.rate,lower=lower,upper=upper,
                        loglik=-optim.out$value)


## 7. Gamma Mixture Distribution 
int.fxn.gamma <- function(theta,x,shape,scale){
  pdf.under.alt.symbolic(x,theta)*dgamma(theta,shape=shape,scale=scale)
}
like.fxn.one.gamma <- function(x,shape,scale){
  integrate(int.fxn.gamma,0,Inf,x=x,shape=shape,scale=scale)$value
}
nll.fxn.gamma <- function(pars,dat.train){
  mean=pars[["mean"]]
  scale=pars[["scale"]]
  shape=mean/scale
  total.log.like <- 0
  for(ii in 1:length(dat.train)){
    tmp.like <- like.fxn.one.gamma(dat.train[[ii]],
                                   shape=shape,scale=scale)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

par0 <- c(mean=1,scale=1)
nll.fxn.gamma(pars=par0,dat.train=dat.train$p)

optim.out <- optim(par=par0,
                   nll.fxn.gamma,
                   dat.train=dat.train$p,
                   method="L-BFGS-B",
                   lower=c(2,0.5),upper=c(4,1))
est.par <- optim.out$par
hess <- hessian(nll.fxn.gamma,optim.out$par,dat.train=dat.train$p)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- est.par - z.crit*sigma
upper <- est.par + z.crit*sigma

gamma.fit <- list(pars=est.par,lower=lower,upper=upper)
gamma.fit$loglik <- -optim.out$value

## 8. Gamma Mixture Distribution, Linear Time 
int.fxn.gamma.lt <- function(theta,x,time.diff,theta0,slope,scale){
  shape <- (theta0+slope*time.diff)/scale
  pdf.under.alt.symbolic(x,theta)*dgamma(theta,shape=shape,scale=scale)
}
like.fxn.one.gamma.lt <- function(x,time.diff,theta0,slope,scale=scale){
  integrate(int.fxn.gamma.lt,0,Inf,x=x,
            time.diff=time.diff,
            theta0=theta0,
            slope=slope,scale=scale)$value
}
nll.fxn.gamma.lt <- function(pars,dat.train){
  theta0=pars[["theta0"]]
  slope=pars[["slope"]]
  scale=pars[["scale"]]
  total.log.like <- 0
  for(ii in 1:nrow(dat.train)){
    tmp.like <- like.fxn.one.gamma.lt(dat.train$p[ii],dat.train$time.diff[ii],
                                    theta0=theta0,slope=slope,scale=scale)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

t0 <- as.numeric(gamma.fit$par["mean"])
par0 <- c(theta0=t0,slope=0.03,scale=gamma.fit$pars[["scale"]])
nll.fxn.gamma.lt(pars=par0,dat.train=dat.train)

optim.out <- optim(par=par0,
                   nll.fxn.gamma.lt,
                   dat.train=dat.train,
                   method="L-BFGS-B",
                   lower=c(2,0.02,0.5),upper=c(4,0.04,1))
est.par <- optim.out$par
hess <- hessian(nll.fxn.gamma.lt,optim.out$par,dat.train=dat.train)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- est.par - z.crit*sigma
upper <- est.par + z.crit*sigma

gamma.fit.lt <- list(pars=est.par,lower=lower,upper=upper,
                        loglik=-optim.out$value)


## Save Data
save(one.theta.nt,one.theta.lt,two.theta.nt,two.theta.lt,
     exp.fit,exp.fit.rate,exp.fit.mean.lt,
     gamma.fit,gamma.fit.lt,
     file="_R_Data_Files/model_fits.RData")

save.image(file="_R_Data_Files/analysis.RData")
