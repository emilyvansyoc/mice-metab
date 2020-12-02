## Final Figure 6 Plot
# Emily Bean, 12/2020

# load packages and get greyscale colors
require(tidyverse)
require(ggpubr)
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

## ---- FIG. 6: Tumor differences between treatments (PA+ER and SED+AL) ----

# run: R/aqueousANOVAs.R lines 9-13, 25-74

# get only PA+ER vs SED+AL
paer <- sigs %>% filter(contrast == "SED+AL-PA+ER")

# filter data
plotdat <- abdat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(!treatmentID == "SED+AL") %>% 
  mutate(treatmentID = factor(treatmentID, ordered = TRUE, levels = c("SED+AL", "PA+AL", "SED+ER", "PA+ER")))


## ---- Plot function ----



fig6 <- function(MetaboliteName, binwidth = 0.2) {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "treatmentID", y = "ab.change", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "treatmentID", width = 0.4, binwidth = binwidth),
                 # change axis titles
                 ylab = "\u0394 SED+AL",
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#969696", "#000000")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              ggtheme = theme_pubr()) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))+
    rremove("xlab") 
  
  return(p1)
  
}

## ---- build Plots ----

# get list of unique metabolites
mets <- unique(plotdat$Metabolite)

# build plots individually
p1 <- fig6(mets[1])
p2 <- fig6(mets[2])
p3 <- fig6(mets[3])
p4 <- fig6(mets[4])
p5 <- fig6(mets[5], binwidth = 0.4) # not sure what's funky with Quin acid but quick fix for now
p6 <- fig6(mets[6])

## arrange and export as one plot
## we want quinolinic acid (met5) to be LAST
g1 <- ggarrange(p1, p2, p3, p4, p6, p5,
                ncol = 2, nrow = 3,
                labels = c("A", "B", "C", "D", "E", "F", "G"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig6.jpeg", plot = g1,  device = "jpeg", dpi = 600, height = 6, width = 8, units = "in")
