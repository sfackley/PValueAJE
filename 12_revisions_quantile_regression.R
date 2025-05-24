rm(list=ls())
library(quantreg)

load(file="_R_Data_Files/full_dataset.RData")

dat$CI.width <- as.numeric(dat$CI_Upper) - as.numeric(dat$CI_Lower)
dat$CI.width <- abs(dat$CI.width)

# Fit the model
qr.model <- rq(CI.width ~ time.diff, data = dat, tau = 0.5)

# Get summary with bootstrapped standard errors
boot.qr <- summary(qr.model, se = "boot", R = 1000)

coefs <- boot.qr$coefficients
ci <- cbind(
  Estimate = coefs[,1],
  Lower_95 = coefs[,1] - qnorm(0.975)*coefs[,2],
  Upper_95 = coefs[,1] + qnorm(0.975)*coefs[,2]
)
ci <- signif(ci,2)

write.csv(ci,file="_Outputs/quantile_reg.csv")

