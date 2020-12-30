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

# get data
dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  # hyphenate instead of _ for treatment ID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# get the average of each SED+AL metabolite

sedal <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # get only SED_AL
  filter(treatmentID == "SED+AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area),
            seSedAl = se(area))

# get change compared to SED+AL
abdat <- dat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  left_join(sedal, by = "Metabolite") %>% 
  # calculate change
  mutate(ab.change = area - avgSedAl) %>% 
  # SED+AL is now 0
  mutate(ab.change = case_when(
    treatmentID %in% "SED+AL" ~ 0,
    treatmentID != "SED+AL" ~ ab.change
  ))%>% 
  select(id, treatmentID, Metabolite, ab.change) 

# perform ANOVA between treatment groups
metabs <- unique(abdat$Metabolite)
pvals <- data.frame()

for(i in 1:length(metabs)) {
  
  # define model
  mod <- aov(ab.change ~ treatmentID, data = filter(abdat, Metabolite == metabs[i]))
  
  # do Tukey posthoc
  tukey <- TukeyHSD(mod)
  
  # get p vals
  ps <- data.frame(Metabolite = metabs[i],
                   contrast = rownames(tukey[1]$treatmentID),
                   pval = round(tukey[1]$treatmentID[,4], 3),
                   Fstat = round(summary(mod)[[1]][1,4], 2),
                   numDF = summary(mod)[[1]]$Df[1],
                   denDF = summary(mod)[[1]]$Df[2],
                   row.names = NULL)
  
  # concatenate
  pvals <- rbind(pvals, ps)
  
}

# get significant 
sigs <- pvals %>% 
  filter(pval < 0.05)


## ---- FIG. 6: Tumor differences between treatments (PA+ER and SED+AL) ----


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
p7 <- fig6(mets[7])

## arrange and export as one plot
## we want quinolinic acid (met5) to be LAST
g1 <- ggarrange(p1, p2, p3, p4, p7, p6, p5,
                ncol = 2, nrow = 4,
                labels = c("A", "B", "C", "D", "E", "F", "G"))

# EXPORT
ggsave(filename = "./data/plots/manuscript-plots/fig6.jpeg", plot = g1,  device = "jpeg", dpi = 600, height = 8, width = 8, units = "in")
