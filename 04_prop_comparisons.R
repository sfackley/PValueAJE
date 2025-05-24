rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_Data_Files/analysis.RData")
source("functions.R")

# Check 
btwn.1and5(dat$p)+over.5(dat$p)+lt.1(dat$p)

# Using the cdf of p-value distributions
# Calculate the factions over time for each model
# Notation: 
# .35 = between 0.03 and 0.05
# .15 = between 0.01 and 0.05
# .5 =  below 0.05 
# .1 = below 0.01
# Derive Time Trends in Fractions & Plot
# Model 2: Basic P Curve with Linear Trend
times <- seq(0,24,0.01)
thetas <- one.theta.lt$pars[["slope"]]*times+
               one.theta.lt$pars[["theta0"]]
mod2.35 <- cdf.under.alt(0.05,thetas)-cdf.under.alt(0.03,thetas)
mod2.15 <- cdf.under.alt(0.05,thetas)-cdf.under.alt(0.01,thetas)
mod2.5 <- cdf.under.alt(0.05,thetas)
mod2.1 <- cdf.under.alt(0.01,thetas)

# Model 4: Two dirac delta with Linear Trend
thetas <- two.theta.lt$pars[["slope"]]*times+
  two.theta.lt$pars[["theta0"]]
null.weight <- two.theta.lt$pars[["null.weight"]]
mod4.35 <- two.hypth.cdf(0.05,thetas,null.weight)-two.hypth.cdf(0.03,thetas,null.weight)
mod4.15 <- two.hypth.cdf(0.05,thetas,null.weight)-two.hypth.cdf(0.01,thetas,null.weight)
mod4.5 <- two.hypth.cdf(0.05,thetas,null.weight)
mod4.1 <- two.hypth.cdf(0.01,thetas,null.weight)

# Model 6: Exponential Mixture with Linear Trend
mean <- exp.fit.mean.lt$pars[["slope"]]*times+
  exp.fit.mean.lt$pars[["theta0"]]
rate <- 1/mean
cdf.exp <- Vectorize(cdf.exp)
mod6.35 <- cdf.exp(0.05,rate)-cdf.exp(0.03,rate)
mod6.15 <- cdf.exp(0.05,rate)-cdf.exp(0.01,rate)
mod6.5 <- cdf.exp(0.05,rate)
mod6.1 <- cdf.exp(0.01,rate)

# Model 8: Gamma Mixture with Linear Trend
thetas <- gamma.fit.lt$pars[["slope"]]*times+
  gamma.fit.lt$pars[["theta0"]]
scale <- gamma.fit.lt$pars[["scale"]]
shape <- thetas/scale
cdf.gamma <- Vectorize(cdf.gamma)
mod8.35 <- cdf.gamma(0.05,shape,scale)-cdf.gamma(0.03,shape,scale)
mod8.15 <- cdf.gamma(0.05,shape,scale)-cdf.gamma(0.01,shape,scale)
mod8.5 <- cdf.gamma(0.05,shape,scale)
mod8.1 <- cdf.gamma(0.01,shape,scale)

# Make a single object of these model-derived quantities to plot
time.trends <- data.frame(times,
                      mod2.35,mod2.15,mod2.1,mod2.5,
                      mod4.35,mod4.15,mod4.1,mod4.5,
                      mod6.35,mod6.15,mod6.1,mod6.5,
                      mod8.35,mod8.15,mod8.1,mod8.5)
time.trends.long <- gather(time.trends, key = "model", value = "Fraction", -times)
time.trends.long %<>% separate_wider_delim(model,".",names=c("Model","Comparison"))
time.trends.long$Comparison <- factor(time.trends.long$Comparison,
                                      levels=c("1","5","15","35"),
                                      labels=c("At or below 0.01","At or below 0.05",
                                               "Between 0.01 and 0.05",
                                               "Between 0.03 and 0.05"))
time.trends.long$Model <- factor(time.trends.long$Model,
                                      levels=c("mod2","mod4","mod6","mod8"),
                                      labels=c("Model 2", "Model 4", "Model 6","Model 8"))

# Epi Journal data to plot
# Now we also plot the raw data by year 
# This is the held out data
# First calcuate the fraction in each of the groupings:
# 0.03-0.05, 0.01-0.05, less than 0.01, less than 0.05
empirical.frxn <- dat %>% filter(save.half) %>% ungroup() %>%
  group_by(Year) %>% 
  summarize(f.35=btwn.3and5(p),
            f.15=btwn.1and5(p),
            f.1=lt.1(p),
            f.5=lt.5(p),
            n=n())
# Calculate the numbers in each of the same groups
empirical.frxn %<>% mutate(n.35=f.35*n,
                     n.15=f.15*n,
                     n.1=f.1*n,
                     n.5=f.5*n) 

# Calculate binomial CIs 
empirical.frxn %<>% rowwise %>% mutate(
  lower.35=binom.test(n.35, n)$conf.int[1],
  lower.15=binom.test(n.15, n)$conf.int[1],
  lower.1=binom.test(n.1, n)$conf.int[1],
  lower.5=binom.test(n.5, n)$conf.int[1]
)
empirical.frxn %<>% mutate(
  upper.35=binom.test(n.35, n)$conf.int[2],
  upper.15=binom.test(n.15, n)$conf.int[2],
  upper.1=binom.test(n.1, n)$conf.int[2],
  upper.5=binom.test(n.5, n)$conf.int[2]
)

# Change data format for easy plotting
empirical.frxn %<>% select(-n)
empirical.frxn %<>% mutate_all(as.numeric)
empirical.long <- gather(empirical.frxn, key = "label", value = "value", -Year)
empirical.long %<>% separate_wider_delim(label,".",names=c("label","Comparison"))
empirical.long %<>% filter(label %in% c("f","lower","upper"))
empirical.to.plot <- spread(empirical.long, key = label, value = value)
empirical.to.plot %<>% mutate(Year=Year+0.5,Model="AJE")

empirical.to.plot$Comparison <- factor(empirical.to.plot$Comparison,
                                      levels=c("1","5","15","35"),
                                      labels=c("At or below 0.01","At or below 0.05",
                                               "Between 0.01 and 0.05",
                                               "Between 0.03 and 0.05"))

# Plotting
colors <- c("black",scales::hue_pal()(4))
time.trends1 <- ggplot(time.trends.long,aes(times+2000,Fraction,group=Model,color=Model))+
  geom_errorbar(data=empirical.to.plot,aes(Year,f,ymin=lower,ymax=upper),width=0)+
  geom_line()+
  geom_line(data=empirical.to.plot,aes(Year,f))+
  facet_wrap(. ~ Comparison,scales = "free")+
  xlab("Year")+ theme_classic()+
  labs(color = NULL)+
  scale_color_manual(name=NULL,values = colors,
                     labels = c(
                       "Data",
                       "Model 2:\nP-Curve", 
                       "Model 4:\nTwo Dirac \u03B4 Mixture", 
                       "Model 6:\nExponential Mixture", 
                       "Model 8:\nGamma Mixture"
                     ),
                     guide = guide_legend(
                       nrow = 2,  # Arrange legend labels into two rows
                       byrow = TRUE  # Ensure the labels fill rows first
                     ))+
  theme(legend.position="bottom", 
        legend.direction="horizontal",
        legend.text = element_text(hjust = 0.5),  # Center-align text within entries
        )

# Save
ggsave("_Outputs/time_trends.jpg",time.trends1,width=6,height=4)
save(time.trends1,file="_R_Data_Files/time_trends1.RData")
