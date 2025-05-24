rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_Data_Files/analysis.RData")
source("functions.R")

## Gamma Mixture Distribution, Linear Time 
int.fxn.gamma.t2 <- function(theta,x,time.diff,theta0,slope,slope.sq,scale){
  shape <- (theta0+slope*time.diff+slope.sq*time.diff^2)/scale
  pdf.under.alt.symbolic(x,theta)*dgamma(theta,shape=shape,scale=scale)
}
like.fxn.one.gamma.t2 <- function(x,time.diff,theta0,slope,slope.sq,scale=scale){
  integrate(int.fxn.gamma.t2,0,Inf,x=x,
            time.diff=time.diff,
            theta0=theta0,
            slope=slope,slope.sq=slope.sq,scale=scale)$value
}
nll.fxn.gamma.sq <- function(pars,dat.train){
  theta0=pars[["theta0"]]
  slope=pars[["slope"]]
  slope.sq=pars[["slope.sq"]]
  scale=pars[["scale"]]
  total.log.like <- 0
  for(ii in 1:nrow(dat.train)){
    tmp.like <- like.fxn.one.gamma.t2(dat.train$p[ii],dat.train$time.diff[ii],
                                      theta0=theta0,
                                      slope=slope,
                                      slope.sq=slope.sq,
                                        scale=scale)
    total.log.like <- total.log.like + log(tmp.like)
  }
  -total.log.like
}

t0 <- as.numeric(gamma.fit.lt$pars["theta0"])
sl <- as.numeric(gamma.fit.lt$pars["slope"])
sc <- as.numeric(gamma.fit.lt$pars["scale"])
par0 <- c(theta0=t0,slope=sl,slope.sq=0,scale=sc)
nll.fxn.gamma.sq(pars=par0,dat.train=dat.train)

optim.out <- optim(par=par0,
                   nll.fxn.gamma.sq,
                   dat.train=dat.train,
                   method="L-BFGS-B",
                   lower=c(2,0.02,-0.01,0.5),upper=c(4,0.04,0.01,1))
loglike.sq.model <- -optim.out$value
loglike.lin.model <- gamma.fit.lt$loglik

# LR p-value for a squared term
p.val.squared.term <- pchisq(2*loglike.sq.model-2*loglike.lin.model,df=1,lower.tail = F)

write(p.val.squared.term,file="_Outputs/pval_squared_term.txt")
