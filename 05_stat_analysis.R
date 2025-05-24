rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_data_files/full_dataset.RData")
source("functions.R")
exp.round <- function(x)round(exp(x),3)

num_by_year <- dat %>% group_by(Year) %>% summarize(n=n())
print(num_by_year,n=100)
# Number of CIs
num.cis <- sum(num_by_year$n)

# Number of abstracts
num.abstracts <- dat$PMID %>% unique %>% length
write("Summary",file="_Outputs/num_cis_abstracts.txt")
write(num.cis,file="_Outputs/num_cis_abstracts.txt",append=T)
write(num.abstracts,file="_Outputs/num_cis_abstracts.txt",append=T)
write("",file="_Outputs/num_cis_abstracts.txt",append=T)

# Calculate Fractions
dat %<>% rowwise() %>%
  mutate(lt.1=lt.1(p),
         lt.5=lt.5(p),
         btwn.1and5=btwn.1and5(p),
         btwn.3and5=btwn.3and5(p))

dat %<>%
  mutate(Journal = as.factor(Journal)) %>%
  mutate(Journal = relevel(as.factor(Journal), ref = "Epi"))


# Models overall adjusted for journal
mod0 <- betareg(p ~ time.diff.scaled+Journal, data = dat,link="logit")
#mod0i <- betareg(p ~ time.diff.scaled+Journal, data = dat,link="logit")
#mod0j <- betareg(p ~ time.diff.scaled*Journal, data = dat,link="logit")
#lrtest(mod0,mod0i)
#lrtest(mod0i,mod0j)
mod1 <- glm(data=dat,lt.1~time.diff.scaled+Journal,family='binomial')
mod5 <- glm(data=dat,lt.5~time.diff.scaled+Journal,family='binomial')
mod15 <- glm(data=dat,btwn.1and5~time.diff.scaled+Journal,family='binomial')
mod35 <-glm(data=dat,btwn.3and5~time.diff.scaled+Journal,family='binomial')

filler <- rep("",3)
s1 <- get.lincom(mod0,model="overall")
s2 <- get.lincom(mod1,model="overall")
s3 <- get.lincom(mod5,model="overall")
s4 <- get.lincom(mod15,model="overall")
s5 <- get.lincom(mod35,model="overall")
reg.coef <- rbind(s1,s2,s3,s4,s5)

# Models by journal
mod0 <- betareg(p ~ time.diff.scaled*Journal, data = dat,link="logit")
mod1 <- glm(data=dat,lt.1~time.diff.scaled*Journal,family='binomial')
mod5 <- glm(data=dat,lt.5~time.diff.scaled*Journal,family='binomial')
mod15 <- glm(data=dat,btwn.1and5~time.diff.scaled*Journal,family='binomial')
mod35 <-glm(data=dat,btwn.3and5~time.diff.scaled*Journal,family='binomial')

s1 <- get.lincom(mod0,method="intercept")
s2 <- get.lincom(mod0)
s3 <- get.lincom(mod1)
s4 <- get.lincom(mod5)
s5 <- get.lincom(mod15)
s6 <- get.lincom(mod35)

reg.coef <- rbind(reg.coef,filler,
      s1,filler,
      s2,filler,
      s3,filler,
      s4,filler,
      s5,filler,
      s6,filler)


# Mixed effects overall
mod0me <- glmmTMB(p ~ time.diff.scaled + (1 | PMID),
                 data = dat %>% filter(p<1&p>0),
                 family = beta_family(link = 'logit'))

mod1me <- glmmTMB(data=dat,lt.1~time.diff.scaled+ (1 | PMID),family='binomial') 
mod5me <- glmmTMB(data=dat,lt.5~time.diff.scaled+ (1 | PMID),family='binomial') 
mod15me <- glmmTMB(data=dat,btwn.1and5~time.diff.scaled+ (1 | PMID),family='binomial')
mod35me <- glmmTMB(data=dat,btwn.3and5~time.diff.scaled+ (1 | PMID),family='binomial')

s1 <- (mod0me %>% confint())[2,c(3,1,2)] %>% exp.round
s2 <- (mod1me %>% confint())[2,c(3,1,2)] %>% exp.round
s3 <- (mod5me %>% confint())[2,c(3,1,2)] %>% exp.round
s4 <- (mod15me %>% confint())[2,c(3,1,2)]%>% exp.round
s5 <- (mod35me %>% confint())[2,c(3,1,2)]%>% exp.round
reg.coef.me <- rbind(s1,s2,s3,s4,s5)

mod0me <- glmmTMB(p ~ time.diff.scaled*Journal + (1 | PMID),
                  data = dat %>% filter(p<1&p>0),
                  family = beta_family(link = 'logit'))

mod1me <- glmmTMB(data=dat,lt.1~time.diff.scaled*Journal+ (1 | PMID),family='binomial') 
mod5me <- glmmTMB(data=dat,lt.5~time.diff.scaled*Journal+ (1 | PMID),family='binomial') 
mod15me <- glmmTMB(data=dat,btwn.1and5~time.diff.scaled*Journal+ (1 | PMID),family='binomial')
mod35me <- glmmTMB(data=dat,btwn.3and5~time.diff.scaled*Journal+ (1 | PMID),family='binomial')

s1 <- get.lincom.me(mod0me,method="intercept")
s2 <- get.lincom.me(mod0me)
s3 <- get.lincom.me(mod1me)
s4 <- get.lincom.me(mod5me)
s5 <- get.lincom.me(mod15me)
s6 <- get.lincom.me(mod35me)

reg.coef.me <- rbind(reg.coef.me,filler,
                  s1,filler,
                  s2,filler,
                  s3,filler,
                  s4,filler,
                  s5,filler,
                  s6,filler)
reg.coef %<>% mutate(`Fixed-Effects Models`=
                       paste(Estimate,", 95% CI: (",
                             `2.5 %`,", ",
                             `97.5 %`,")",sep=""))
reg.coef %<>% select(`Fixed-Effects Models`)
reg.coef.me %<>% mutate(`Mixed-Effects Models`=
                       paste(Estimate,", 95% CI: (",
                             `2.5 %`,", ",
                             `97.5 %`,")",sep=""))
reg.coef.me %<>% select(`Mixed-Effects Models`)

reg.coef <- cbind(reg.coef,reg.coef.me)
row.names <- rownames(reg.coef)
numeric.row.names <- grepl("^\\d+$", row.names)
reg.coef[numeric.row.names,] <- ""

write.csv(reg.coef,"_Outputs/stat_coefs.csv")
write.csv(
  dat %>% group_by(Journal) %>% summarize(N=n()),
  file="_Outputs/ns_by_journal.csv")

frac.by.period <- dat %>% group_by(Year>2012) %>%
  summarize(frac.03=signif(btwn.3and5(p)*100,2)) %>%
  print(n=50) 
write.table(frac.by.period,file="_Outputs/frac_by_period.txt")

save.image(file="_R_data_files/stat_analysis.RData")
