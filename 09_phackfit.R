rm(list=ls())
load(file="_R_Data_Files/analysis.RData")
load(file="_R_Data_Files/stat_analysis.RData")
source('functions.R')


do.sim <- function(mean,scale,shape,times,mean.attempts.bl,slope,n.ps){
    p.vals <- c() 
    for(ii in 1:n.ps){
      mean.attempts <- mean.attempts.bl + slope*times[ii]
      theta <- rgamma(1,shape=shape,scale=scale)
      p <- get.random(theta)
      max.attempts <- rpois(1,mean.attempts) 
      attempts <- 0
      while(p>0.05 & attempts<max.attempts){
        p <- get.random(theta)
        attempts <- attempts+1
      }
      p.vals <- c(p.vals,p) 
    }
  p.vals
}

do.many.sim <- function(n.sims,mean,scale,shape,times,mean.attempts.bl,slope,n.ps){
  out <- matrix(ncol=5,nrow=0)
  for(ii in 1:n.sims){
    if(ii%%10==0)print(ii)
    p.vals <- do.sim(mean,scale,shape,times,mean.attempts.bl,slope,n.ps)
    dat.sim <- dat
    p.vals[p.vals<1e-10] <- 1e-10
    p.vals[p.vals>1-1e-10] <- 1e-10
    dat.sim$p <- p.vals
    mod <- betareg(p ~ time.diff.scaled, data = dat.sim,link="logit")
    time.coef <- mod$coefficients$mean[["time.diff.scaled"]]
    tmp.out <- c("betareg"=time.coef) 
    dat.sim %<>% rowwise() %>%
      mutate(lt.1=lt.1(p),
             lt.5=lt.5(p),
             btwn.1and5=btwn.1and5(p),
             btwn.3and5=btwn.3and5(p))
    mod1 <- glm(data=dat.sim,lt.1~time.diff.scaled,family='binomial')
    mod5 <- glm(data=dat.sim,lt.5~time.diff.scaled,family='binomial')
    mod15 <- glm(data=dat.sim,btwn.1and5~time.diff.scaled,family='binomial')
    mod35 <-glm(data=dat.sim,btwn.3and5~time.diff.scaled,family='binomial')
    tmp.out["lr1"] <- mod1$coefficients[["time.diff.scaled"]]
    tmp.out["lr5"] <- mod5$coefficients[["time.diff.scaled"]]
    tmp.out["lr15"] <- mod15$coefficients[["time.diff.scaled"]]
    tmp.out["lr35"] <- mod35$coefficients[["time.diff.scaled"]]
    out <- rbind(out,tmp.out)
  }
  out
}

get.sim.summary <- function(sim.slopes){
  sim.slopes %<>% as.data.frame %>% mutate_if(is.numeric,exp)
  sim.summary <- data.frame(
    mean=apply(sim.slopes,2,mean),
    se=apply(sim.slopes,2,function(x)sd(x)/sqrt(length(x))),
    min=apply(sim.slopes,2,min),
    max=apply(sim.slopes,2,max)
  )
  sim.summary %<>% mutate(lower=mean-qnorm(0.975)*se)
  sim.summary %<>% mutate(upper=mean+qnorm(0.975)*se)
  sim.summary %<>% mutate_if(is.numeric,function(x)round(x,3))
  
  sim.summary %<>% mutate(Mean=paste(mean,", 95% CI: (",lower,", ",upper,")",sep=""))
  sim.summary %<>% mutate(Range=paste("[",min,", ",max,"]",sep=""))
  sim.summary %<>% select(Mean,Range)
  sim.summary
}

# No Time Trend
n.ps <- nrow(dat)
n.sims <- 100
mean <- gamma.fit$pars["mean"]
scale <- gamma.fit$pars["scale"]
shape <- mean/scale
times <- dat$time.diff
mean.attempts.bl <- 5
slope <- 0  

set.seed(12345)
sim.slopes.nt <- do.many.sim(n.sims,mean,scale,shape,times,mean.attempts.bl,slope,n.ps)

# With Time Trend -- decreased p-hacking
mean.attempts.bl <- 5
slope <- - 0.4

set.seed(12345)
sim.slopes.lt <- do.many.sim(n.sims,mean,scale,shape,times,mean.attempts.bl,slope,n.ps)

# With Time Trend -- increased p-hacking
mean.attempts.bl <- 5
slope <- 0.4 

set.seed(12345)
sim.slopes.inc.ph <- do.many.sim(n.sims,mean,scale,shape,times,mean.attempts.bl,slope,n.ps)


tab1 <- get.sim.summary(sim.slopes.nt)
tab2 <- get.sim.summary(sim.slopes.lt)
tab3 <- get.sim.summary(sim.slopes.inc.ph)

sim.phack.table <- rbind(rep("",ncol(tab1)),tab1,
                         rep("",ncol(tab1)),tab2,
                         rep("",ncol(tab1)),tab3)

write.csv(sim.phack.table,file="_Outputs/sim_phack_table.csv")
save.image(file="_R_Data_Files/phack.RData")
