rm(list=ls())
select <- dplyr::select

load(file="_R_Data_Files/model_fits.RData")

# Make a dataframe with all model fits:
# Model 1
m1 <- as.data.frame(one.theta.nt)
m1.text <- data.frame(Model="1. One d",Parameters="d")
m1 <- cbind(m1.text,m1)
m1$n.pars <- nrow(m1)
# Model 2
m2 <- as.data.frame(one.theta.lt)
m2.text <- data.frame(Model="2. One d, linear trend",
                      Parameters=c("Time trend","d at 1/1/2000"))
m2 <- cbind(m2.text,m2)
m2 <- m2[2:1,]
m2$n.pars <- nrow(m2)

# Model 3
m3 <- as.data.frame(two.theta.nt)
m3.text <- data.frame(Model="3. Two d",Parameters=c("Mean d","Null weight"))
m3 <- cbind(m3.text,m3)
m3$n.pars <- nrow(m3)

# Model 4
m4 <- as.data.frame(two.theta.lt)
m4.text <- data.frame(Model="4. Two d, linear trend",
                      Parameters=c("Time trend", "Mean d at 1/1/2000","Null weight"))
m4 <- cbind(m4.text,m4)
m4 <- m4[c(2,3,1),]
m4$n.pars <- nrow(m4)

# Model 5
m5 <- as.data.frame(exp.fit)
m5.text <- data.frame(Model="5. Exponential mixture",
                      Parameters=c("Mean d"))
m5 <- cbind(m5.text,m5)
m5$n.pars <- nrow(m5)

# Model 6
m6 <- as.data.frame(exp.fit.mean.lt)
m6.text <- data.frame(Model="6. Exponential mixture, linear trend",
                      Parameters=c("Mean d at 1/1/2000","Time trend"))
m6 <- cbind(m6.text,m6)
m6$n.pars <- nrow(m6)

# Model 7
m7 <- as.data.frame(gamma.fit)
m7.text <- data.frame(Model="7. Gamma mixture",
                      Parameters=c("Mean d","Scale"))
m7 <- cbind(m7.text,m7)
m7$n.pars <- nrow(m7)

# Model 8
m8 <- as.data.frame(gamma.fit.lt)
m8.text <- data.frame(Model="8. Gamma mixture, linear trend",
                      Parameters=c("Mean d at 1/1/2000","Time trend","Scale"))
m8 <- cbind(m8.text,m8)
m8$n.pars <- nrow(m8)

# Combine:
full.table <- rbind(m1,m2,m3,m4,m5,m6,m7,m8) 
power.fun <- function(theta){
  z.crit <- qnorm(0.975)
  power <- 1+pnorm(-z.crit-theta)-pnorm(z.crit-theta)
  round(power*100,1)
}
full.table %<>% mutate(power=power.fun(pars))
full.table %<>% mutate(fp1=", 95% CI: (")
full.table %<>% mutate(LowerP=power.fun(lower))
full.table %<>% mutate(fp2=", ")
full.table %<>% mutate(UpperP=power.fun(upper))
full.table %<>% mutate(fp3=")")

# Include 95% CIs in the table
full.table$Est <- sub("^0+", "", signif(full.table$pars,3)) 
full.table$f1 <- ", 95% CI: ("
full.table$Lower <- sub("^0+", "", signif(full.table$lower,3)) 
full.table$f2 <- ", "
full.table$Upper <- sub("^0+", "", signif(full.table$upper,3)) 
full.table$f3 <- ")"

# Pull out LL's for later: 
like.tab <- full.table %>% select(Model,loglik,n.pars) %>% distinct()
rownames(like.tab) <- NULL

# Drop unneeded columns 
full.table %<>% unite(Estimate,Est:f3,sep="")
full.table %<>% unite(Power,power:fp3,sep="")
full.table %<>% select(-pars,-lower,-upper,-loglik,-n.pars)
rownames(full.table) <- NULL
full.table %<>% select(-Power, everything())
full.table %<>% mutate(Power=ifelse(Power>70,Power,""))
full.table %<>% rename(`Power (%)`=Power)
# Write table of model parameters to CSV
write.csv(full.table,file="_Outputs/mod_params_table3.csv",row.names=F)

# Round values 
# Include the log likelihoods in the table
# Make column names: model 1, model 2, etc.
like.tab %<>% mutate(AIC=2*n.pars-2*loglik)
like.tab %<>% mutate(loglik=round(loglik,1))
like.tab %<>% mutate(AIC=round(AIC,1))
like.tab %<>% rename(`Log Likelihood`=loglik,
                     `Number of Parameters`=n.pars)

# Save as csv
write.csv(like.tab,file="_Outputs/likelihood_table_S3.csv",row.names=F)
