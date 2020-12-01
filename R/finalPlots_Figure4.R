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

## ---- FIG 4: Plasma at Day 35 ----

# run: R/aqueousANOVAs.R line 9-13, 101-173

# also get marginally significant metabolites (p < 0.1)

# get significant and marginally significant
msig <- pvals %>% filter(contrast == "SED+AL-PA+ER") %>% filter(pval < 0.1)

# filter data
plotdat <- d35 %>% 
  semi_join(msig, by = "Metabolite") %>% 
  filter(!treatmentID == "SED+AL") %>% 
  mutate(treatmentID = factor(treatmentID, ordered = TRUE, levels = c("PA+AL", "SED+ER", "PA+ER")))

## ---- plot ----

# define plot
p <- ggbarplot(data = plotdat, x = "treatmentID", y = "ab.change", fill = "#D9D9D9",
               #facet.by = "Metabolite", scales = "free",
               add = c("mean_se", "dotplot"), add.params = list(fill = "treatmentID", width = 0.4, binwidth = 0.2),
               # change axis titles
               xlab = "treatmentID", ylab = "\u0394 SED+AL",
               # facet wrap
               facet.by = "Metabolite") +
  scale_fill_manual(values = treatmentGreys) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 

g <- ggpar(p, legend = "none",
            xlab = FALSE,
            ggtheme = theme_pubr()) +
  # make facet wrap title background white
  theme(strip.background = element_rect(
    fill="white", linetype=0
  ), 
  strip.text.x = element_text(size =  12))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig4.jpeg", plot = g,  device = "jpeg", dpi = 600, height = 4, width = 7, units = "in")
