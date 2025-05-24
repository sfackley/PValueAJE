rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select
recode <- dplyr::recode

# Read in data set produced in Python
dat <- read.delim("../../extracted_data.tsv", header = TRUE, sep = "\t")

# Re-code names to abbreviations
journal_recode <- c(
  "Am J Epidemiol" = "AJE",
  "Epidemiology" = "Epi",
  "Eur J Epidemiol" = "EJE",
  "Int J Epidemiol" = "IJE"
)
dat %<>% mutate(Journal = recode(Journal, !!!journal_recode))

# Create a ratio column that indicates whether something is a ratio measure
dat$Ratio <- (dat$Measure_Category=="ratio")

# Coerce everything to numeric
dat %<>% mutate(
  Estimate = as.numeric(Point_Estimate),
  Lower=as.numeric(CI_Lower),
  Upper=as.numeric(CI_Upper)
)

# Logarithmize estimates and CIs if on the ratio scale
dat %<>% mutate(Estimate=ifelse(Ratio,
                                log(Estimate),
                                Estimate))
dat %<>% mutate(Lower=ifelse(Ratio,
                                log(Lower),
                                Lower))
dat %<>% mutate(Upper=ifelse(Ratio,
                                log(Upper),
                                Upper))

dropped.due.to.log <- 
(is.nan(dat$Estimate)|is.nan(dat$Lower)|is.nan(dat$Upper)) %>% sum()
write(dropped.due.to.log,file="_Outputs/log_dropped.txt")

dat %<>% mutate(Calc.Arg=abs(Estimate)/(abs(Upper-Lower)/2/qnorm(0.975)))
dat %<>% mutate(logp=pnorm(Calc.Arg,log.p=T,lower.tail=F)+log(2))
dat %<>% mutate(Calc.P.Val=2-2*pnorm(Calc.Arg))
dat %<>% mutate(p=exp(logp))

# Check to make sure exp(logp) comes out to the p-value
any((dat$p-dat$Calc.P.Val)>0.0000001,na.rm=T)

# Data cleaning so p-values are suitable for beta-regression
dat %<>% mutate(p=ifelse(p<1e-10,1e-10,p)) # Reassign p-values close to zero
dat %<>% mutate(p=ifelse(p>1-1e-10,1-1e-10,p)) # Reassign p-values close to 1
dat %<>% drop_na(p)

# Fix dates
# If date is missing append the 15th of the month (mid-month)
dat %<>% mutate(Date=ifelse(!grepl("\\d$",Date),paste0(Date,"15"),Date))
# If month is also missing, assume July 2nd (midyear) 
dat %<>% mutate(Date = ifelse(!grepl("[A-Za-z]", Date), paste0(Date, "July 2"), Date))
dat %<>% mutate(Date=as.Date(Date,format = "%Y %b %d"))
# New variable is time.diff and gives time from the midpoint of the period
dat %<>% mutate(time.diff=as.numeric(interval(as.Date("2012-07-01"),Date)/dyears(1)))
# Rescale time diff for per decade changes
dat %<>% mutate(time.diff.scaled=time.diff/10)

# Remove any dates before 2000 or after 2024
dat <- dat %>%
  filter(Date > as.Date("2000-01-01") & Date < as.Date("2024-12-31"))
dat$Year <- floor(year(dat$Date))

# Save dataset for next script
save(dat,file="_R_Data_Files/full_dataset.RData")
epi <- dat %>% filter(Journal=="Epi")
save(epi,file="_R_Data_Files/epi_dataset.RData")
aje <- dat %>% filter(Journal=="AJE")
save(aje,file="_R_Data_Files/aje_dataset.RData")
eje <- dat %>% filter(Journal=="EJE")
save(eje,file="_R_Data_Files/eje_dataset.RData")
ije <- dat %>% filter(Journal=="IJE")
save(ije,file="_R_Data_Files/ije_dataset.RData")

