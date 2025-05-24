load(file="_R_Data_Files/phack.RData")
load(file="_R_Data_Files/time_trends1.RData")

n.ps <- nrow(dat)
n.sims <- 100
mean <- gamma.fit$pars["mean"]
scale <- gamma.fit$pars["scale"]
shape <- mean/scale
times <- dat$time.diff
mean.attempts.bl <- 1
slope <- 0  

#
mean.attempts.bl <- 1
model.ps <- 
  do.sim(mean,scale,shape,times,mean.attempts.bl,slope,n.ps)

#
mean.attempts.bl <- 5
p.hacking <- 
  do.sim(mean,scale,shape,times,mean.attempts.bl,slope,n.ps)

#
mean.attempts.bl <- 1
more.power <- 
  do.sim(mean,scale*2,shape,times,mean.attempts.bl,slope,n.ps)

# Combine the vectors into a data frame
data <- data.frame(
  value = c(model.ps, p.hacking, more.power),
  source = factor(
    c(rep("Model Parameters", length(model.ps)),
      rep("P Hacking", length(p.hacking)),
      rep("Increased Power", length(more.power)))
  )
)
data$source <- fct_relevel(data$source, "Model Parameters", 
                           "P Hacking", 
                           "Increased Power")

colors <- c(scales::hue_pal()(4))
library(ggplot2)
library(cowplot)

# Main plot
first_plot <- ggplot(data, aes(x = value, fill = source)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.05,  # Binwidth for the main plot
                 alpha = 0.7, 
                 position = "identity", 
                 color = "black") +  # Outlined bars
  scale_fill_manual(name = NULL, values = colors) +
  labs(
    x = "Values",
    y = "Density",
    fill = "Source"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  facet_wrap(~ source, scales = "free_x")

# Adjust binwidth only for the second plot
second_plot <- ggplot(data, aes(x = value, fill = source)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.01,  # Smaller binwidth for more bars
                 alpha = 0.7, 
                 position = "identity", 
                 color = "black",breaks=seq(0,1,0.01)) +  # Outlined bars
  scale_fill_manual(name = NULL, values = colors) +
  labs(
    x = "P Value",
    y = "Square-Root Density",
    fill = "Source"
  ) +
  scale_y_sqrt() +  # Logarithmic scale for y-axis
  theme_classic() +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    panel.spacing = unit(1, "lines")  # Increase space between plots
  ) +
  facet_wrap(~ source, scales = "free_x") +
  coord_cartesian(xlim = c(0.0, 0.0765))  # Zoomed range


# Combine main plot and inset

combined <- time.trends1 / second_plot +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(ncol = 1, heights = c(0.45, 0.25))

ggsave("_Outputs/figure3.jpg",combined,width=6,height=8,)
ggsave("_Outputs/figure3.tiff",combined,width=6,height=8,)

