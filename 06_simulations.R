rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_data_files/analysis.RData")
n.sims <- 100
set.seed(2000)
source("functions.R")
dat.test <- dat[save.half,]
n.ps <- nrow(dat.test)

# Model 1: One Theta, No time
theta <- one.theta.nt$pars
mod1.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%10==0)print(ii)
  mod1.sim[,ii] <- get.n.random(n.ps,theta) 
}
summarize.sim.dat(mod1.sim)

# Model 2: One Theta, Linear time
thetas <- one.theta.lt$pars[["slope"]]*dat.test$time.diff+
  one.theta.lt$pars[["theta0"]]
mod2.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%10==0)print(ii)
  mod2.sim[,ii] <- get.random(thetas) 
}
summarize.sim.dat(mod2.sim)

# Model 3: Two Thetas, No time
theta <- two.theta.nt$pars[["theta"]]
null.weight <- two.theta.nt$pars[["null.weight"]]
mod3.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%10==0)print(ii)
  mod3.sim[,ii] <- get.n.random.2hypoth(n.ps,theta,null.weight) 
}
summarize.sim.dat(mod3.sim)

# Model 4: Two Thetas, Linear time
thetas <- two.theta.lt$pars[["slope"]]*dat.test$time.diff+
  two.theta.lt$pars[["theta0"]]
null.weight <- two.theta.lt$pars[["null.weight"]]
mod4.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%10==0)print(ii)
  mod4.sim[,ii] <- get.random.2hypoth(thetas,null.weight) 
}
summarize.sim.dat(mod4.sim)

# Model 5: Exponential, No time
rate <- 1/exp.fit$pars[["theta"]]
mod5.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%5==0)print(ii)
  mod5.sim[,ii] <- get.n.random.exp(n.ps,rate) 
}
summarize.sim.dat(mod5.sim)

# Model 6: Exponential, Linear time
theta <- exp.fit.mean.lt$pars[["slope"]]*dat.test$time.diff+
  exp.fit.mean.lt$pars[["theta0"]]
mod6.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%5==0)print(ii)
  mod6.sim[,ii] <- get.random.exp(rate=1/theta) 
}
summarize.sim.dat(mod6.sim)  

# Model 7: Gamma, No time
mean <- gamma.fit$pars[["mean"]]
scale <- gamma.fit$pars[["scale"]]
shape <- mean/scale
mod7.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%5==0)print(ii)
  mod7.sim[,ii] <- get.n.random.gamma(n.ps,shape,scale) 
}
summarize.sim.dat(mod7.sim)

# Model 8: Gamma, Linear time
theta <- gamma.fit.lt$pars[["slope"]]*dat.test$time.diff+
  gamma.fit.lt$pars[["theta0"]]
scale <- gamma.fit.lt$pars[["scale"]]
shape <- theta/scale
mod8.sim <- matrix(nrow=n.ps,ncol=n.sims)
for(ii in 1:n.sims){
  if(ii%%5==0)print(ii)
  mod8.sim[,ii] <- get.random.gamma(shape,scale) 
}
summarize.sim.dat(mod8.sim) 

# Save results
f1 <- matrix(c("","","",""),1)
f2 <- matrix(c("","","",""),1)
f3 <- matrix(c("","","",""),1)
f4 <- matrix(c("","","",""),1)
f5 <- matrix(c("","","",""),1)
f6 <- matrix(c("","","",""),1)
f7 <- matrix(c("","","",""),1)
f8 <- matrix(c("","","",""),1)

rownames(f1) <- "Model 1"
rownames(f2) <- "Model 2"
rownames(f3) <- "Model 3"
rownames(f4) <- "Model 4"
rownames(f5) <- "Model 5"
rownames(f6) <- "Model 6"
rownames(f7) <- "Model 7"
rownames(f8) <- "Model 8"
colnames(f1) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f2) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f3) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f4) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f5) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f6) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f7) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
colnames(f8) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")

x <- rbind(f1,
           summarize.sim.dat(mod1.sim),
           f2,
           summarize.sim.dat(mod2.sim),
           f3,
           summarize.sim.dat(mod3.sim),
           f4,
           summarize.sim.dat(mod4.sim),
           f5,
           summarize.sim.dat(mod5.sim),
           f6,
           summarize.sim.dat(mod6.sim),
           f7,
           summarize.sim.dat(mod7.sim),
           f8,
           summarize.sim.dat(mod8.sim))
write.csv(x,file="_Outputs/pmodel_summary.csv")

med.time <- median(dat.test$time.diff)+12.5
x <- rbind(f2,
           summarize.sim.dat.by.time(mod2.sim,med.time),
           f4,
           summarize.sim.dat.by.time(mod4.sim,med.time),
           f6,
           summarize.sim.dat.by.time(mod6.sim,med.time),
           f8,
           summarize.sim.dat.by.time(mod8.sim,med.time))
write.csv(x,file="_Outputs/pmodel_summary_median.csv")

x <- rbind(f2,
           summarize.sim.dat.by.time(mod2.sim,dichot.time=16),
           f4,
           summarize.sim.dat.by.time(mod4.sim,dichot.time=16),
           f6,
           summarize.sim.dat.by.time(mod6.sim,dichot.time=16),
           f8,
           summarize.sim.dat.by.time(mod8.sim,dichot.time=16))
write.csv(x,file="_Outputs/pmodel_summary_16.csv")

x <- rbind(f2,
           summarize.sim.dat.by.time(mod2.sim,dichot.time=12.5),
           f4,
           summarize.sim.dat.by.time(mod4.sim,dichot.time=12.5),
           f6,
           summarize.sim.dat.by.time(mod6.sim,dichot.time=12.5),
           f8,
           summarize.sim.dat.by.time(mod8.sim,dichot.time=12.5))
write.csv(x,file="_Outputs/pmodel_summary_12.csv")

save.image(file="_R_Data_Files/sims.RData")
