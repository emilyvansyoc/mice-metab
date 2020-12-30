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

# get data
dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  # hyphenate instead of _ for treatment ID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-")) %>% 
  # get only plasma
  filter(str_detect(id, "plasma")) %>% 
  # remove day 21
  filter(!str_detect(id, "_plasmaD21")) %>% 
  # add Time column
  mutate(Time = case_when(
    str_detect(id, "plasmaD7") ~ "Day 7",
    str_detect(id, "plasmaD35") ~ "Day 35"
  ))

## Important: Lactic acid is not detected in plasma; need to remove 
dat <- dat %>% filter(!Metabolite == "L-Lactic acid")

# perform t test between Day 7 and Day 35
# get unique metabolites and treatments
trt <- unique(dat$treatmentID)
mets <- unique(dat$Metabolite)

pvals <- data.frame()

for(j in 1:length(trt)) {
  
  for(i in 1:length(mets)) {
    
    # define model
    mod <- t.test(area ~ Time, data = filter(dat, treatmentID == trt[j] & Metabolite == mets[i]))
    
    # get p vals
    ps <- data.frame(Metabolite = mets[i],
                     treatmentID = trt[j],
                     pval = round(mod$p.value, 3),
                     tstat = round(mod$statistic, 3),
                     row.names = NULL)
    
    # concatenate
    pvals <- rbind(pvals, ps)
  }
}

# get significant 
sigs <- pvals %>% 
  drop_na() %>% 
  filter(pval < 0.05)

# get only PA+ER
paer <- sigs %>% filter(treatmentID == "PA+ER")


# get table of mean +/- SE for results section
df <- dat %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  group_by(treatmentID, Metabolite, Time) %>% 
  summarize(avgchange = round(mean(area), 3),
            sechange = round(se(area), 3)) %>% 
  pivot_wider(names_from = c(treatmentID, Time), values_from = c(avgchange, sechange))

# write to table
#write.table(df, "./data/plasma-ttest-day7vsday35.txt", sep = "\t", row.names = FALSE)

# which metabolites are significant in each time
t <- sigs %>% select(-tstat) %>% pivot_wider(names_from = treatmentID, values_from = pval)

## ---- FIG 3: Plasma over time ----

# filter data
plotdat <- dat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(treatmentID == "PA+ER") %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 35")))


## ---- Plot function ----

# since there are 14 metabolites, build a function to make them individually


fig3 <- function(MetaboliteName, ylab = "Relative concentration") {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "Time", y = "area", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "Time", width = 0.4, binwidth = 0.2),
                 # change axis titles
                 xlab = "Time", ylab = ylab,
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#525252")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              ggtheme = theme_pubr(), ylab = ylab) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))+
    rremove("x.text") +
    rremove("xlab") 
  
  return(p1)
  
}

# define a similar function to include labels for the bottom plots
fig3Labs <- function(MetaboliteName, ylab = "Relative concentration") {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "Time", y = "area", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "Time", width = 0.4, binwidth = 0.2),
                 # change axis titles
                 xlab = "Time", ylab = ylab,
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#525252")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              ggtheme = theme_pubr(), ylab = ylab) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
  
}

  
## ---- make plots ----

mets <- unique(plotdat$Metabolite)

p1 <- fig3(mets[1])
p2 <- fig3(mets[2], ylab = FALSE)
p3 <- fig3(mets[3])
p4 <- fig3(mets[4], ylab = FALSE)
p5 <- fig3(mets[5])
p6 <- fig3(mets[6], ylab = FALSE)
p7 <- fig3(mets[7])
p8 <- fig3(mets[8], ylab = FALSE)
p9 <- fig3(mets[9])
p10 <- fig3(mets[10], ylab = FALSE)
p11 <- fig3(mets[11])
p12 <- fig3(mets[12], ylab = FALSE)
# the last 2 need labels since they will be at the bottom
p13 <- fig3Labs(mets[13])
p14 <- fig3Labs(mets[14], ylab = FALSE)



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
ggsave(filename = "./data/plots/manuscript-plots/fig3-part4.jpeg", plot = g4,  device = "jpeg", dpi = 600, height = 4.5, width = 5, units = "in")
