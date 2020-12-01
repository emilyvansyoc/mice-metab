## Final Figure 3 Plot
# Emily Bean, 12/2020

# load packages and get greyscale colors
require(tidyverse)
require(ggpubr)
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

## ---- FIG 3: Plasma over time ----

# run: ./R/aqueousANOVAs.R lines 9-13, 101-135, 230-266

# filter data
plotdat <- tdat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(treatmentID == "PA+ER") %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 35")))


## ---- Plot function ----

# since there are 14 metabolites, build a function to make them individually


fig3 <- function(MetaboliteName) {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "Time", y = "ab.change", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "Time", width = 0.4, binwidth = 0.2),
                 # change axis titles
                 xlab = "Time", ylab = "\u0394 SED+AL",
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#525252")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              xlab = FALSE,
              ggtheme = theme_pubr()) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
  
  }

  
## ---- make plots ----

mets <- unique(plotdat$Metabolite)

p1 <- fig3(mets[1])
p2 <- fig3(mets[2])
p3 <- fig3(mets[3])
p4 <- fig3(mets[4])
p5 <- fig3(mets[5])
p6 <- fig3(mets[6])
p7 <- fig3(mets[7])
p8 <- fig3(mets[8])
p9 <- fig3(mets[9])
p10 <- fig3(mets[10])
p11 <- fig3(mets[11])
p12 <- fig3(mets[12])
p13 <- fig3(mets[13])
p14 <- fig3(mets[14])



## add all plots together and export 
# exporting all in one makes the dots too small - break it up
g1 <- ggarrange(p1, 
                p2, 
                p3,
                p4,
                ncol = 2, nrow = 2,
                labels = c("A", "B", "C", "D"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig3-part1.jpeg", plot = g1,  device = "jpeg", dpi = 600, height = 4, width = 5, units = "in")


##2
g2 <- ggarrange(p5, 
                p6, 
                p7,
                p8,
                ncol = 2, nrow = 2,
                labels = c("E", "F", "G", "H"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig3-part2.jpeg", plot = g2,  device = "jpeg", dpi = 600, height = 4, width = 5, units = "in")

##3
g3 <- ggarrange(p9, 
                p10, 
                p11,
                p12,
                ncol = 2, nrow = 2,
                labels = c("I", "J", "K", "L"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig3-part3.jpeg", plot = g3,  device = "jpeg", dpi = 600, height = 4, width = 5, units = "in")

##4
g4 <- ggarrange(p13, 
                p14, 
                ncol = 2, nrow = 2,
                labels = c("M", "N"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig3-part4.jpeg", plot = g4,  device = "jpeg", dpi = 600, height = 4, width = 5, units = "in")
