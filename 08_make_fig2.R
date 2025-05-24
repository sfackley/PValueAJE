rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_Data_Files/stat_analysis.RData")
load(file="_R_Data_Files/full_dataset.RData")
load(file="_R_Data_Files/analysis.RData")


# Table 1
pref.order <- c("Epi", "AJE", "EJE", "IJE")
dat$Journal <- factor(dat$Journal, levels = pref.order)
dat %<>% mutate(Year=floor(time.diff+2012))
summary <- dat %>% group_by(Journal) %>%
  summarize(
    `N of Abstracts`=length(unique(PMID)),
    `N of CIs`=n(),
    `Mean of P`=signif(mean(p),3),
    med=signif(median(logp,na.rm=T),3),
    low.iqr=signif(quantile(logp,0.25),3),
    high.iqr=signif(quantile(logp,0.75),3)
  )
overall.summary <- 
  dat %>% 
  summarize(
    `N of Abstracts`=length(unique(PMID)),
    `N of CIs`=n(),
    `Mean of P`=signif(mean(p),3),
    med=signif(median(logp,na.rm=T),3),
    low.iqr=signif(quantile(logp,0.25),3),
    high.iqr=signif(quantile(logp,0.75),3)
  )
overall.summary$Journal="Overall"
summary <- rbind(summary,overall.summary)
summary %<>% mutate(`Median of log P (IQR)`=paste(med,", (",low.iqr,", ", high.iqr,")",sep=""))
summary %<>% select(-med,-low.iqr,-high.iqr)

write.csv(summary,file="_Outputs/table1.csv")

dat$Year <- floor(dat$time.diff+2012.5)
summary.supp <- dat %>% group_by(Journal,Year) %>%
  summarize(
    `N of Abstracts`=length(unique(PMID)),
    `N of CIs`=n(),
    `Mean of P`=signif(mean(p),3),
    med=signif(median(logp,na.rm=T),3),
    low.iqr=signif(quantile(logp,0.25),3),
    high.iqr=signif(quantile(logp,0.75),3)
  )
summary.supp %<>% mutate(`Median of log P (IQR)`=paste(med,", (",low.iqr,", ", high.iqr,")",sep=""))
summary.supp %<>% select(-med,-low.iqr,-high.iqr)
write.csv(summary.supp,file="_Outputs/tableS1.csv")

summary.to.plot <- dat %>% group_by(Journal,Year) %>%
  summarize(
    `N of Abstracts`=length(unique(PMID)),
    `N of CIs`=n(),
    `Mean of P`=signif(mean(p),3),
    med=signif(median(p,na.rm=T),3),
    low.iqr=signif(quantile(p,0.25),3),
    high.iqr=signif(quantile(p,0.75),3)
  )

overall <- 
  dat %>% group_by(Year) %>%
  summarize(
    Journal="Overall",
    `N of Abstracts`=length(unique(PMID)),
    `N of CIs`=n(),
    `Mean of P`=signif(mean(p)),
    med=(median(p,na.rm=T)),
    low.iqr=(quantile(p,0.25)),
    high.iqr=(quantile(p,0.75))) %>%
  select("Year","Journal",everything())

summary.to.plot <- rbind(summary.to.plot,overall)

# By journal by year plot
set.seed(12345)
pref.order <- c("Epi", "AJE", "EJE", "IJE","Overall")
summary.to.plot$Journal <- factor(summary.to.plot$Journal, levels = pref.order)
colors <- c(scales::hue_pal()(4),"black")
lwds <- c(rep(0.5, 4), 1)
by.journal.plot <- ggplot(summary.to.plot,
                          aes(Year,med,color=Journal,
                              group=Journal,
                              linewidth=Journal,
                              ymin=low.iqr,
                              ymax=high.iqr))+
  geom_line()+
  geom_errorbar(aes(ymin=low.iqr,ymax=high.iqr),
                alpha=0.5,width=0,
                position=position_jitter(width=0.2)) +
  scale_color_manual(name=NULL,values = colors,
                     labels = c("Epi","AJE", "EJE", "IJE", "Overall"),
                     guide = guide_legend(
                       override.aes = list(
                         linewidth = lwds)))+
  theme_minimal()+
  ylab(expression(paste("Median ", italic("P"), " (IQR), Square-Root Scale")))+
  theme(legend.position = c(0.49, 0.75), legend.direction="horizontal")+
  scale_linewidth_manual(values = lwds,
                         labels = c("Epi","AJE", "EJE", "IJE", "Overall"),
                         guide=NULL)+
  scale_y_sqrt(breaks=seq(0,0.4,0.05)^2,minor_breaks=NULL)
by.journal.plot

# Histogram
main.hist <- ggplot(dat,aes(p)) +
  geom_histogram(color='darkgrey',fill='white',breaks=seq(0,1,0.01),alpha=0.1) +
  scale_y_sqrt(breaks=seq(0,80,10)^2,minor_breaks=NULL)+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  theme_minimal()+xlab(expression(paste(italic("P")," Value")))+ylab("Count, Square-Root Scale")+
  geom_vline(xintercept=0.05,col='black')+
  geom_vline(xintercept=0.01,col='black')

inset.hist <- ggplot(dat %>% filter(p<=0.1), aes(p)) +
  geom_histogram(color='darkgrey',fill='white',breaks=seq(0,.1,0.001),alpha=0.1) +
  scale_x_continuous(breaks=c(0.01,0.05,0.03,0.02, 0.04, 0.07))+
  coord_cartesian(xlim=c(0.008,0.055),ylim=c(0,350))+
  scale_y_continuous(minor_breaks=NULL)+
  labs(title = "Zoomed Inset") +
  theme_minimal()+xlab("")+ylab("Count, Linear Scale")+
  theme(plot.background = element_rect(fill = "white", color = "black", linewidth = 1))+
  geom_vline(xintercept=0.05,col='black')+
  geom_vline(xintercept=0.01,col='black')  

# Combine the plots 
#main.hist + inset.hist + plot_layout(guides = "collect")
histv1 <- main.hist + 
  annotation_custom(ggplotGrob(inset.hist), xmin = 0.2, xmax = 0.9, ymin = 12, ymax = 95)
ggsave("_Outputs/histogram.jpg",histv1,width=6,height=4)


theta <- one.theta.nt$pars
fit1 <- function(x,t=theta)pdf.under.alt.symbolic(x,t)

theta <- two.theta.nt$pars[["theta"]]
null.weight <- two.theta.nt$pars[["null.weight"]]
fit3 <- function(x,t=theta,nw=null.weight)two.hypth(x,t,nw)

rate <- 1/exp.fit$pars[["theta"]]
fit5 <- function(x,r=rate)like.fxn.one.exp(x,r)
fit5 <- Vectorize(fit5)

mean <- gamma.fit$pars[["mean"]]
scale <- gamma.fit$pars[["scale"]]
shape <- mean/scale
fit7 <- function(x,sh=shape,sc=scale)like.fxn.one.gamma(x,sh,sc)
fit7 <- Vectorize(fit7)

# Histogram
colors <- scales::hue_pal()(4)
lwd <- 0.5
main.hist <- ggplot(dat.train,aes(p)) +
  stat_function(fun = fit1, aes(color = colors[1]),linewidth = lwd,n=10000,show.legend=TRUE)+
  stat_function(fun = fit3, aes(color = colors[2]), linewidth = lwd,n=10000,show.legend=TRUE)+
  stat_function(fun = fit5, aes(color = colors[3]), linewidth = lwd,n=10000,show.legend=TRUE)+
  stat_function(fun = fit7, aes(color = colors[4]), linewidth = lwd,n=10000,show.legend=TRUE)+
  geom_histogram(aes(y = ..density..),
                 color='darkgrey',fill='white',breaks=seq(0,1,0.01),alpha=0.1) +
  scale_y_sqrt(breaks=seq(1,10)^2,minor_breaks=NULL,lim=c(0,70))+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  theme_minimal()+xlab(expression(paste(italic("P")," Value")))+ylab("Density, Square-Root Scale")+
  scale_color_manual(name=NULL,values = colors,
                     labels = c(
                       "Model 2:\nP-Curve", 
                       "Model 4:\nTwo Dirac \u03B4 Mixture", 
                       "Model 6:\nExponential Mixture", 
                       "Model 8:\nGamma Mixture"
                     ),
                     guide = "legend")+
  theme(legend.position=c(0.535, 0.3), 
        legend.direction="horizontal",
        legend.text = element_text(size = 7))+
  geom_vline(xintercept=0.05,col='black')+
  geom_vline(xintercept=0.01,col='black')
main.hist

inset.hist <- ggplot(dat.train, aes(p))+
  geom_histogram(aes(y = ..density..),
                 color='darkgrey',fill='white',breaks=seq(0,.1,0.001),alpha=0.5) +
  scale_x_continuous(breaks=c(0.01,0.05,0.03,0.02, 0.04, 0.07))+                                              # Zoom without cutting values
  coord_cartesian(xlim=c(0.008,0.055))+
  scale_y_continuous(minor_breaks=NULL,lim=c(0,20))+
  labs(title = "Zoomed Inset") +
  theme_minimal()+xlab("")+ylab("Density, Linear Scale")+
  theme(plot.background = element_rect(fill = "white", color = "black", linewidth = 1))+
  stat_function(fun = fit1, color = colors[1], linewidth = lwd,n=10000)+
  stat_function(fun = fit3, color = colors[2], linewidth = lwd,n=10000)+
  stat_function(fun = fit5, color = colors[3], linewidth = lwd,n=10000)+
  stat_function(fun = fit7, color = colors[4], linewidth = lwd,n=10000)+
  geom_vline(xintercept=0.05,col='black')+
  geom_vline(xintercept=0.01,col='black')
inset.hist
# Combine the plots 
#main.hist + inset.hist + plot_layout(guides = "collect")
histv2 <- main.hist + 
  annotation_custom(ggplotGrob(inset.hist), xmin = 0.2, xmax = 0.9, ymin = 3, ymax = 9)
histv2
ggsave("_Outputs/histogram_v2.jpg",histv2,width=6,height=4)

combined <- histv1 / histv2 / by.journal.plot +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(ncol = 1, heights = c(0.45, 0.45,0.45))

ggsave("_Outputs/histogram_combined.jpg",combined,width=6,height=10,)
ggsave("_Outputs/figure2.tiff",combined,width=6,height=10,)
ggsave("_Outputs/figure2.eps",combined,width=6,height=10,)

