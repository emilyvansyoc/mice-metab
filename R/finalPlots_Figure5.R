## Final Figure 4 Plot
# Emily Bean, 12/2020

# load packages and get greyscale colors
require(tidyverse)
require(ggpubr)
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

## ---- FIG. 5: Tumor corr to tumor volume ----

# run: R/aqueousTumorVol.R lines 9-66

# make plot function to plot individually

fig4 <- function(MetaboliteName) {
  
  # define data
  dat <- plotdat %>% filter(Metabolite == MetaboliteName)
  
  # create plot 
  plot <- ggscatter(dat, x = "area", y = "cm3",
                  # add regression line
                  add = "reg.line", conf.int = FALSE,
                  repel = TRUE,
                  # add Pearons R on the plot
                  cor.coef = TRUE,
                  cor.method = "pearson",
                  add.params = list(color = "black",
                                    fill = "black")) +
    facet_wrap(~ Metabolite, scales = "free") +
    #stat_cor(method = "pearson") +
    labs(x = "Relative concentration", y = expression(paste("Tumor volume ", cm^3)))  +
    theme(strip.background = element_rect(
      fill="white", linetype=0
    ))
  
  return(plot)
}

## ---- Positive slopes ----

# get significant metabs
pos <- sigs %>% filter(sign == "pos" & r2 > 0.2)

# filter data
plotdat <- df %>% 
  semi_join(pos, by = "Metabolite")

# get list of metabolites
metpos <- unique(plotdat$Metabolite)

## plot individually
p1 <- fig4(metpos[1])
p2 <- fig4(metpos[2])
p3 <- fig4(metpos[3])
p4 <- fig4(metpos[4])
p5 <- fig4(metpos[5])
p6 <- fig4(metpos[6])


## ---- Negative Slopes ----

# get significant metabs
neg <- sigs %>% filter(sign == "neg" & r2 > 0.2)

# filter data
plotdat <- df %>% 
  semi_join(neg, by = "Metabolite")

# get list of metabolites
metneg <- unique(plotdat$Metabolite)

# plot individually
p7 <- fig4(metneg[1])
p8 <- fig4(metneg[2])
p9 <- fig4(metneg[3])

## arrange and export as one plot
g1 <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,
                ncol = 3, nrow = 3, 
                labels = c("A", "B", "C", "D", "E", "F", "G", "H", "J"))
# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig5.jpeg", plot = g1,  device = "jpeg", dpi = 600, height = 7, width = 8.5, units = "in")

