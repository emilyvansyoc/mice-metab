### Finalized Figure 2 Plot
# Emily Bean, 11/2020


# load packages and get greyscale colors
require(tidyverse)
require(ggpubr)
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))


### ---- Figure 2: Plasma regression to tumor volume ----

# run: ./R/aqueousTumorVol.R lines 9-35, lines 74-97 (chunks "read data" and "plasma regressions")

# plot positive and negative slopes separately
# and only r2 > 0.2
pos <- sigs %>% filter(r2 > 0.2 & sign == "pos")
neg <- sigs %>% filter(r2 > 0.2 & sign == "neg")

## ---- Positive slopes ----

# create the plotting dataframe
plotdat2pos <- df %>% 
  semi_join(pos, by = "Metabolite")

# get list of unique metabolites to plot individually
mets2a <- unique(plotdat2pos$Metabolite)

# create plot 
p2a <- ggscatter(filter(plotdat2pos, Metabolite == mets2a[1]), x = "area", y = "cm3",
               # add regression line
               add = "reg.line", conf.int = FALSE,
               # shape and fill by Time
               shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

p2b <- ggscatter(filter(plotdat2pos, Metabolite == mets2a[2]), x = "area", y = "cm3",
                # add regression line
                add = "reg.line", conf.int = FALSE,
                # shape and fill by Time
                shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

p2c <- ggscatter(filter(plotdat2pos, Metabolite == mets2a[3]), x = "area", y = "cm3",
                # add regression line
                add = "reg.line", conf.int = FALSE,
                # shape and fill by Time
                shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

## add all plots together and export
g1 <- ggarrange(p2a, p2b, p2c, ncol = 3, nrow = 1, labels = c("A", "B", "C"),
          common.legend = TRUE)
# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/Fig2-positive.jpeg", plot = g1,  device = "jpeg", dpi = 300, height = 4, width = 8, units = "in")



## ----  Negative slopes ----

# create the plotting dataframe
plotdat2neg <- df %>% 
  semi_join(neg, by = "Metabolite")

# get list of unique metabolites to plot individually
mets2b <- unique(plotdat2neg$Metabolite)

# create plot 
p2d <- ggscatter(filter(plotdat2neg, Metabolite == mets2b[1]), x = "area", y = "cm3",
                 # add regression line
                 add = "reg.line", conf.int = FALSE,
                 # shape and fill by Time
                 shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

p2e <- ggscatter(filter(plotdat2neg, Metabolite == mets2b[2]), x = "area", y = "cm3",
                 # add regression line
                 add = "reg.line", conf.int = FALSE,
                 # shape and fill by Time
                 shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

p2f <- ggscatter(filter(plotdat2neg, Metabolite == mets2b[3]), x = "area", y = "cm3",
                 # add regression line
                 add = "reg.line", conf.int = FALSE,
                 # shape and fill by Time
                 shape = "Time", color = "Time", palette = c("black", "darkgrey"),
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

## add all plots together and export
g2 <- ggarrange(p2d, p2e, p2f, ncol = 3, nrow = 1, labels = c("D", "E", "F"),
          font.label = list(size = 14, color = "black", face = "bold"),
          common.legend = TRUE) 
# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/Fig2-negative.jpeg", plot = g2,  device = "jpeg", dpi = 300, height = 4, width = 8, units = "in")


