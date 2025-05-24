cdf.under.alt <- function(x,theta){
  x1 <- x/2
  x2 <- 1-x/2
  pnorm(qnorm(x1)+theta)+pnorm(qnorm(x2)+theta,lower.tail=F)
}

two.hypth.cdf <- function(x,theta,null.weight){
  null.weight*cdf.under.alt(x,0)+
    (1-null.weight)*cdf.under.alt(x,theta)
}

get.random <- function(theta){
  ran.unif <- runif(1)
  find.root <- function(x)cdf.under.alt(x,theta)-ran.unif
  uniroot(find.root,c(0,1))$root
}
get.random <- Vectorize(get.random)

get.n.random <- function(n,theta){
  vec <- c()
  for(ii in 1:n){
    vec <- c(vec,get.random(theta))
  }
  vec
}

get.random.2hypoth <- function(theta,null.weight){
  ran.unif <- runif(1)
  find.root <- function(x)two.hypth.cdf(x,theta,null.weight)-ran.unif
  uniroot(find.root,c(0,1))$root
}
get.random.2hypoth <- Vectorize(get.random.2hypoth)

get.n.random.2hypoth <- function(n,theta,null.weight){
  vec <- c()
  for(ii in 1:n){
    vec <- c(vec,get.random.2hypoth(theta,null.weight))
  }
  vec
}

int.fxn.exp <- function(theta,x,rate){
  cdf.under.alt(x,theta)*dexp(theta,rate=rate)
}
cdf.exp <- function(x,rate){
  integrate(int.fxn.exp,0,Inf,x=x,rate=rate)$value
}

get.random.exp <- function(rate){
  ran.unif <- runif(1)
  find.root <- function(x)cdf.exp(x,rate)-ran.unif
  uniroot(find.root,c(0,1))$root
}
get.random.exp <- Vectorize(get.random.exp)

get.n.random.exp <- function(n,rate){
  vec <- c()
  for(ii in 1:n){
    vec <- c(vec,get.random.exp(rate))
  }
  vec
}

int.fxn.gamma <- function(theta,x,shape,scale){
  cdf.under.alt(x,theta)*dgamma(theta,shape=shape,scale=scale)
}
cdf.gamma <- function(x,shape=shape,scale=scale){
  integrate(int.fxn.gamma,0,Inf,x=x,shape=shape,scale=scale)$value
}

get.random.gamma <- function(shape,scale){
  ran.unif <- runif(1)
  find.root <- function(x)cdf.gamma(x,shape,scale)-ran.unif
  uniroot(find.root,c(0,1))$root
}
get.random.gamma <- Vectorize(get.random.gamma)

get.n.random.gamma <- function(n,shape,scale){
  vec <- c()
  for(ii in 1:n){
    vec <- c(vec,get.random.gamma(shape,scale))
  }
  vec
}

btwn.1and5 <- function(x){
  sum(x<=0.05&x>0.01)/length(x)
}
btwn.3and5 <- function(x){
  sum(x<=0.05&x>=0.03)/length(x)
}
over.5 <- function(x){
  sum(x>0.05)/length(x)
}
lt.1 <- function(x){
  sum(x<=0.01)/length(x)
}

lt.5 <- function(x){
  sum(x<=0.05)/length(x)
}

summarize.sim.dat <- function(sim.dat){
  sim.dat.summary <- data.frame(
    btwn.1and5=apply(sim.dat,2,btwn.1and5),
    over.5=apply(sim.dat,2,over.5), 
    lt.1=apply(sim.dat,2,lt.1),
    btwn.3and5=apply(sim.dat,2,btwn.3and5),
    lt.5=apply(sim.dat,2,lt.5)
  )
  s1 <- sim.dat.summary %>% 
    summarize(mean=mean(btwn.1and5),
              min=min(btwn.1and5),
              max=max(btwn.1and5),
              dat=btwn.1and5(dat.test$p))
  s1p <- sim.dat.summary %>% 
    summarize(mean=mean(btwn.3and5),
              min=min(btwn.3and5),
              max=max(btwn.3and5),
              dat=btwn.3and5(dat.test$p))
  s2 <- sim.dat.summary %>% 
    summarize(mean=mean(over.5),
              min=min(over.5),
              max=max(over.5),
              dat=over.5(dat.test$p))
  s3 <- sim.dat.summary %>% 
    summarize(mean=mean(lt.1),
              min=min(lt.1),
              max=max(lt.1),
              dat=lt.1(dat.test$p))
  s4 <- sim.dat.summary %>% 
    summarize(mean=mean(lt.5),
              min=min(lt.5),
              max=max(lt.5),
              dat=lt.5(dat.test$p))
  tab <- rbind(s1,s1p,s2,s3,s4) %>% round(.,3)
  rownames(tab) <- c("Between 0.01 and 0.05",
                     "Between 0.03 and 0.05",
                     "Over 0.05", 
                     "Less than 0.01",
                     "All less than 0.05 ")
  colnames(tab) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
  tab
}

summarize.sim.dat.by.time <- function(sim.dat,dichot.time=NULL){
  #print(dichot.time)
  dichot.time <- dichot.time - 12.5
  if(is.null(dichot.time)){
    early.ps <- dat.test$p[dat.test$time.diff.scaled<=median(dat.test$time.diff.scaled)]
    late.ps <- dat.test$p[dat.test$time.diff.scaled>median(dat.test$time.diff.scaled)]
    sim.dat.early <- sim.dat[dat.test$time.diff.scaled<=median(dat.test$time.diff.scaled),]
    sim.dat.late <- sim.dat[dat.test$time.diff.scaled>median(dat.test$time.diff.scaled),]
  }else{
    early.ps <- dat.test$p[dat.test$time.diff<=dichot.time]
    late.ps <- dat.test$p[dat.test$time.diff>dichot.time]
    sim.dat.early <- sim.dat[dat.test$time.diff<=dichot.time,]
    sim.dat.late <- sim.dat[dat.test$time.diff>dichot.time,]
  }
  sim.dat.summary.early <- data.frame(
    btwn.1and5=apply(sim.dat.early,2,btwn.1and5),
    over.5=apply(sim.dat.early,2,over.5), 
    lt.1=apply(sim.dat.early,2,lt.1),
    lt.5=apply(sim.dat.early,2,lt.5),
    btwn.3and5=apply(sim.dat.early,2,btwn.3and5)
  )
  sim.dat.summary.late <- data.frame(
    btwn.1and5=apply(sim.dat.late,2,btwn.1and5),
    over.5=apply(sim.dat.late,2,over.5), 
    lt.1=apply(sim.dat.late,2,lt.1),
    lt.5=apply(sim.dat.late,2,lt.5),
    btwn.3and5=apply(sim.dat.late,2,btwn.3and5)
  )
  s1.early <- sim.dat.summary.early %>% 
    summarize(mean=mean(btwn.1and5),
              min=min(btwn.1and5),
              max=max(btwn.1and5),
              dat=btwn.1and5(early.ps))
  s1p.early <- sim.dat.summary.early %>% 
    summarize(mean=mean(btwn.3and5),
              min=min(btwn.3and5),
              max=max(btwn.3and5),
              dat=btwn.3and5(early.ps))
  s2.early <- sim.dat.summary.early %>% 
    summarize(mean=mean(over.5),
              min=min(over.5),
              max=max(over.5),
              dat=over.5(early.ps))
  s3.early <- sim.dat.summary.early %>% 
    summarize(mean=mean(lt.1),
              min=min(lt.1),
              max=max(lt.1),
              dat=lt.1(early.ps))
  s4.early <- sim.dat.summary.early %>% 
    summarize(mean=mean(lt.5),
              min=min(lt.5),
              max=max(lt.5),
              dat=lt.5(early.ps))
  s1.late <- sim.dat.summary.late %>% 
    summarize(mean=mean(btwn.1and5),
              min=min(btwn.1and5),
              max=max(btwn.1and5),
              dat=btwn.1and5(late.ps))
  s1p.late <- sim.dat.summary.late %>% 
    summarize(mean=mean(btwn.3and5),
              min=min(btwn.3and5),
              max=max(btwn.3and5),
              dat=btwn.3and5(late.ps))
  s2.late <- sim.dat.summary.late %>% 
    summarize(mean=mean(over.5),
              min=min(over.5),
              max=max(over.5),
              dat=over.5(late.ps))
  s3.late <- sim.dat.summary.late %>% 
    summarize(mean=mean(lt.1),
              min=min(lt.1),
              max=max(lt.1),
              dat=lt.1(late.ps))
  s4.late <- sim.dat.summary.late %>% 
    summarize(mean=mean(lt.5),
              min=min(lt.5),
              max=max(lt.5),
              dat=lt.5(late.ps))
  tab <- rbind(s1.early,s1.late,
               s1p.early,s1p.late,
               s2.early,s2.late,
               s3.early,s3.late,
               s4.early,s4.late) %>% 
    round(.,3)
  r.names <- c("Between 0.01 and 0.05","Between 0.03 and 0.05",
               "Over 0.05", "Less than 0.01","All less than 0.05")
  rownames(tab) <- c(paste("Early:",r.names),paste("Late:",r.names))[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(tab) <- c("Mean Fraction","Minimum Fraction", "Maximum Fraction","Fraction in Epi Journals")
  tab
}


# Lincom function 
get.lincom <- function(mod,method="slope",model="individual"){
  if(method=="intercept"){
    out <- lincom(mod,c("(Intercept)",
                        "(Intercept)+JournalAJE",
                        "(Intercept)+JournalEJE",
                        "(Intercept)+JournalIJE")) %>%
      as.data.frame() %>%
      mutate_all(as.numeric) %>%
      mutate_if(is.numeric,function(x)round(exp(x),3)) %>%
      select(-ncol(.)) %>% select(-ncol(.)) 
    rownames(out) <- c("Epidemiology","AJE","EJE","IJE")
  }
  if(model=="overall"){
    out <- lincom(mod,c("time.diff.scaled")) %>%
      as.data.frame() %>%
      mutate_all(as.numeric) %>%
      mutate_if(is.numeric,function(x)round(exp(x),3)) %>%
      select(-ncol(.)) %>% select(-ncol(.)) 
    rownames(out) <- c("Overall")
  }
  if(!method=="intercept"&!model=="overall"){
    out <- lincom(mod,c("time.diff.scaled",
                        "time.diff.scaled+time.diff.scaled:JournalAJE",
                        "time.diff.scaled+time.diff.scaled:JournalEJE",
                        "time.diff.scaled+time.diff.scaled:JournalIJE")) %>%
      as.data.frame() %>%
      mutate_all(as.numeric) %>%
      mutate_if(is.numeric,function(x)round(exp(x),3)) %>%
      select(-ncol(.)) %>% select(-ncol(.))
    rownames(out) <- c("Epidemiology","AJE","EJE","IJE")
  }
  out
}

get.lincom.me <- function(mod,method="slope"){
  sum.mod <- summary(mod)$coefficients$cond
  z.crit <- qnorm(0.975)
  if(method=="intercept"){
    out <- data.frame(Estimate=
                        rep(sum.mod["(Intercept)","Estimate"],4))
    out$Estimate[2:4] <- out$Estimate[2:4]+
        sum.mod[c("JournalAJE","JournalEJE","JournalIJE"),"Estimate"]
    out$var <- rep(sum.mod["(Intercept)","Std. Error"]^2)
    out$var[2:4] <- out$var[2:4]+
      sum.mod[c("JournalAJE","JournalEJE","JournalIJE"),"Std. Error"]^2
    out$`2.5 %` <- out$Estimate - z.crit*sqrt(out$var)
    out$`97.5 %` <- out$Estimate + z.crit*sqrt(out$var)
    out %<>% as.data.frame() %>%
      mutate_all(as.numeric) %>%
      mutate_if(is.numeric,function(x)round(exp(x),3)) %>%
      select(c("Estimate","2.5 %","97.5 %"))
  }else{
    out <- data.frame(Estimate=
                        rep(sum.mod["time.diff.scaled","Estimate"],4))
    out$Estimate[2:4] <- out$Estimate[2:4]+
      sum.mod[c("time.diff.scaled:JournalAJE",
                "time.diff.scaled:JournalEJE",
                "time.diff.scaled:JournalIJE"),"Estimate"]
    out$var <- rep(sum.mod["time.diff.scaled","Std. Error"]^2)
    out$var[2:4] <- out$var[2:4]+
      sum.mod[c("time.diff.scaled:JournalAJE",
                "time.diff.scaled:JournalEJE",
                "time.diff.scaled:JournalIJE"),"Std. Error"]^2
    out$`2.5 %` <- out$Estimate - z.crit*sqrt(out$var)
    out$`97.5 %` <- out$Estimate + z.crit*sqrt(out$var)
    out %<>% as.data.frame() %>%
      mutate_all(as.numeric) %>%
      mutate_if(is.numeric,function(x)round(exp(x),3)) %>%
      select(c("Estimate","2.5 %","97.5 %"))
  }
  rownames(out) <- c("Epidemiology","AJE","EJE","IJE")
  out
}

