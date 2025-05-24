rm(list=ls())
mutate <- dplyr::mutate
select <- dplyr::select

load(file="_R_Data_Files/sims.RData")
load(file="_R_Data_Files/analysis.RData")
source("functions.R")

summarize.sim.by.year <- function(mod,fun){
  mod <- as.data.frame(mod)
  mod$Year <- dat.test$Year 
  fraction <- mod %>%
    group_by(Year) %>%
    summarize(across(starts_with("V"), ~ fun(.)))
  result <- fraction %>%
    mutate(
      mean = rowMeans(select(., starts_with("V")), na.rm = TRUE),
      min = apply(select(., starts_with("V")), 1, min, na.rm = TRUE),
      max = apply(select(., starts_with("V")), 1, max, na.rm = TRUE)
    )
  result %>% select(-contains("V"))
}

mods <- c(2,4,6,8)
mod.list <- paste("mod",mods,".sim",sep="")
label.list <- paste("Model",mods)
fun.list <- c("btwn.3and5","btwn.1and5","lt.1","lt.5")
comparison <- c("35","15","1","5")

sim.summaries <- data.frame()
for(mod in mod.list){
  for(fun in fun.list){
    expr.string <- paste("summarize.sim.by.year(",mod,",",fun,")",sep="")
    expr.parsed <- parse(text=expr.string)
    sim.sum <- eval(expr.parsed)
    sim.sum$Model <- label.list[mod.list==mod]
    sim.sum$Comparison <- comparison[fun.list==fun]
    sim.summaries <- rbind(sim.summaries,sim.sum)
  }
}
sim.summaries$Comparison <- factor(sim.summaries$Comparison,
                                 levels=c("1","5","15","35"),
                                 labels=c("At or below 0.01","At or below 0.05",
                                          "Between 0.01 and 0.05",
                                          "Between 0.03 and 0.05"))
sim.summaries %<>% filter(!Year==2024)

# Empirical data to plot
empirical.frxn <- dat.test %>% group_by(Year) %>% #
  summarize(f.35=btwn.3and5(p),
            f.15=btwn.1and5(p),
            f.1=lt.1(p),
            f.5=lt.5(p))
empirical.frxn %<>% filter(!Year==2024)
empirical.long <- gather(empirical.frxn, key = "label", value = "value", -Year)
empirical.long %<>% separate_wider_delim(label,".",names=c("label","Comparison"))
empirical.to.plot <- spread(empirical.long, key = label, value = value)
empirical.to.plot %<>% mutate(Year=Year,Model="AJE")
empirical.to.plot$Comparison <- factor(empirical.to.plot$Comparison,
                                 levels=c("1","5","15","35"),
                                 labels=c("At or below 0.01","At or below 0.05",
                                          "Between 0.01 and 0.05",
                                          "Between 0.03 and 0.05"))

# Plotting
set.seed(12345) # set seed for the jitter
colors <- c("black",scales::hue_pal()(4))
time.trends2 <- ggplot(sim.summaries,aes(Year,mean,group=Model,color=Model))+
  geom_line(lwd=1)+
  geom_errorbar(aes(ymin=min,ymax=max),
                alpha=0.5,width=0,
                position=position_jitter(width=0.2))+
  facet_wrap(. ~ Comparison,scales = "free")+
  theme_classic()+
  geom_line(data=empirical.to.plot,aes(Year,f))+
  xlab("Year")+ 
  labs(color = NULL)+
  ylab("Fraction")+
  scale_color_manual(name=NULL,values = colors,
                     labels = c("AJE","Model 2", "Model 4", "Model 6", "Model 8"),
                     guide = NULL)+
  scale_x_continuous(breaks=seq(2000,2025,5),limits=c(2000,2025))+
  theme(legend.position="bottom", legend.direction="horizontal")
ggsave("_R_Data_Files/time_trends2.jpg",time.trends2,width=6,height=4)

load(file="_R_Data_Files/time_trends1.RData")
combined <- time.trends1 / time.trends2 +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(ncol = 1, heights = c(0.45, 0.45))
ggsave("_Outputs/time_trends_combined.jpg",combined,width=6,height=8,)
