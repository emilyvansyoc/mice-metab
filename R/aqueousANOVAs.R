
## Aqueous metabolite ANOVAs on tumor and plasma tissue

# ---- get data ----
require(tidyverse)

# define standard error function
se <- function(x) sqrt(var(x)/length(x))

dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  # hyphenate instead of _ for treatment ID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# set ggplot theme 
theme_set(theme_bw())

# get colors
source("./R/RColorBrewer.R")

## ---- tumor ----

# get the average of each SED-AL metabolite

sedal <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # get only SED_AL
  filter(treatmentID == "SED-AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area),
            seSedAl = se(area))

# get change compared to SED-AL
abdat <- dat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  left_join(sedal, by = "Metabolite") %>% 
  # calculate change
  mutate(ab.change = area - avgSedAl) %>% 
  # SED-AL is now 0
  mutate(ab.change = case_when(
    treatmentID %in% "SED-AL" ~ 0,
    treatmentID != "SED-AL" ~ ab.change
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
                   pval = tukey[1]$treatmentID[,4],
                   row.names = NULL)
  
  # concatenate
  pvals <- rbind(pvals, ps)
  
}

# get significant 
sigs <- pvals %>% 
  filter(pval < 0.05)

# get table of mean +/- SE for results section
df <- abdat %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  group_by(treatmentID, Metabolite) %>% 
  summarize(avgchange = round(mean(ab.change), 3),
            sechange = round(se(ab.change), 3)) %>% 
  pivot_wider(names_from = treatmentID, values_from = c(avgchange, sechange))

# write to table
#write.table(df, "./data/changefromSEDAL-tumor-ANOVA.txt", sep = "\t", row.names = FALSE)

## PLOTS --- make one plot for each of the metabolites, then put together in PPT

# re-order treatments: EX-AL, SED-ER, EX-ER
plotdat <- abdat %>% 
  filter(!treatmentID == "SED-AL") %>% 
  mutate(treatment = factor(treatmentID, ordered = TRUE, levels = c("EX-AL", "SED-ER", "EX-ER"))) %>% 
  semi_join(sigs, by = "Metabolite")

## do this in a loop
met <- unique(plotdat$Metabolite)

for(i in 1:length(met)) {
  # create plot
  plot <- ggplot(data = filter(plotdat, Metabolite == met[i]), aes(x = treatment, y = ab.change, fill = treatment)) +
    geom_boxplot() +
    geom_jitter(width = 0.08, size = 2) +
    scale_fill_manual(values = treatmentIntcols) +
    theme_classic() +
    labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment") +
    ggtitle(as.character(met[i]))
  
  # save
  ggsave(filename = paste0("./data/plots/indiv-boxplot-tumor-ANOVA-", as.character(met[i]), ".png"), plot = plot, dpi = 600, height = 4.35, width = 7.32, units = "in")
  
  # print progress
  cat(as.character(met[i]))
}

# for manuscript: need number of replicates (mice) in each treatment
reps <- abdat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  group_by(treatmentID) %>% 
  summarize(reps = length(unique(id)))



## PLOTS: Weight effects under exercise (EX-ER vs EX-AL)
wt <- c("5-Thymidylic acid", "Acetylphosphate", "ADP", "D-Glucose", "Hydroxproline", "Indole-3-carboxylic acid",
        "L-Alanine", "L-Arginine", "Succinic acid")
wtdat <- abdat %>% 
  filter(Metabolite %in% wt) %>% 
  filter(!treatmentID %in% "SED-AL")

ggplot(data = wtdat, aes(x = treatmentID, y = ab.change, fill = treatmentID)) +
  geom_boxplot() +
  facet_wrap(~Metabolite) +
  scale_fill_manual(values = treatmentIntcols) +
  labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment")

# save 
#ggsave(filename = "./data/plots/changeSEDAL-tumor-weighteffects.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

## PLOTS: Weight effects under sedentary (SED-ER vs SED-AL)
wt <- c("L-Histidine", "Quinolinic acid")

wtdat <- abdat %>% 
  filter(Metabolite %in% wt) %>% 
  filter(!treatmentID %in% "SED-AL")

ggplot(data = wtdat, aes(x = treatmentID, y = ab.change, fill = treatmentID)) +
  geom_boxplot() +
  facet_wrap(~Metabolite) +
  scale_fill_manual(values = treatmentIntcols) +
  labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment")

## save 
#ggsave(filename = "./data/plots/changeSEDAL-tumor-weighteffectsSEDENTARY.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

## PLOTS: EX-ER vs SED-AL
ex <- c("L-Cystathionine", "L-Glutamine", "L-Lysine", "Quinolinic acid", "Taurine", "D-Glucose", "L-Histadine")

exdat <- abdat %>% 
  filter(Metabolite %in% ex) %>% 
  filter(!treatmentID %in% "SED-AL") %>% 
  group_by(treatmentID, Metabolite) %>% 
  summarize(mean = mean(ab.change),
            se = se(ab.change))

ggplot(data = filter(exdat, treatmentID == "EX-ER"), aes(x = Metabolite, y = mean, fill = treatmentID)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2)) +
  scale_fill_manual(values = treatmentIntcols) +
  labs(x = "Metabolite", y = "Change from SED-AL", fill = "Treatment")

## save 
#ggsave(filename = "./data/plots/changeSEDAL-tumor-EX-ER.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

## ---- plasma ----

# get only plasma and calculate Time column
plasma <- dat %>% 
  filter(str_detect(id, "_plasma")) %>% 
  # make Time column
  mutate(Time = case_when(
    str_detect(id, "_plasmaD7") ~ "Day 7",
    str_detect(id, "_plasmaD21") ~ "Day 21",
    str_detect(id, "_plasmaD35") ~ "Day 35"
  ))

# get the average of each SED-AL metabolite

sedal <- plasma %>% 
  # get only SED_AL
  filter(treatmentID == "SED-AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area),
            seSedAl = se(area))

# get change compared to SED-AL
abdat <- plasma %>% 
  left_join(sedal, by = "Metabolite") %>% 
  # calculate change
  mutate(ab.change = area - avgSedAl) %>% 
  # SED-AL is now 0
  mutate(ab.change = case_when(
    treatmentID %in% "SED-AL" ~ 0,
    treatmentID != "SED-AL" ~ ab.change
  ))%>% 
  select(id, treatmentID, Time, Metabolite, ab.change) 

### 1. ANOVA WITHIN D35 (are there metabolic changes at the end of tumor progression?)

# get only d35
d35 <- abdat %>% 
  filter(Time == "Day 35") %>% 
  mutate(treatmentID = factor(treatmentID))

# do model loop
metabs <- unique(d35$Metabolite)
pvals <- data.frame()

for(i in 1:length(metabs)) {
  
  # define model
  mod <- aov(ab.change ~ treatmentID, data = filter(d35, Metabolite == metabs[i]))
  
  # do Tukey posthoc
  tukey <- TukeyHSD(mod)
  
  # get p vals
  ps <- data.frame(Metabolite = metabs[i],
                   contrast = rownames(tukey[1]$treatmentID),
                   pval = tukey[1]$treatmentID[,4],
                   row.names = NULL)
  
  # concatenate
  pvals <- rbind(pvals, ps)
  
}

# get significant 
sigs <- pvals %>% 
  filter(pval < 0.05)

# get table of mean +/- SE for results section
df <- d35 %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  group_by(treatmentID, Metabolite) %>% 
  summarize(avgchange = round(mean(ab.change), 3),
            sechange = round(se(ab.change), 3)) %>% 
  pivot_wider(names_from = treatmentID, values_from = c(avgchange, sechange))

# write to table
#write.table(df, "./data/changefromSEDAL-plasma-d35-ANOVA.txt", sep = "\t", row.names = FALSE)

### PLOTS --- make one plot for each of the metabolites, then put together in PPT

# re-order treatments: EX-AL, SED-ER, EX-ER
plotdat <- d35 %>% 
  filter(!treatmentID == "SED-AL") %>% 
  mutate(treatment = factor(treatmentID, ordered = TRUE, levels = c("EX-AL", "SED-ER", "EX-ER"))) %>% 
  semi_join(sigs, by = "Metabolite")

## do this in a loop
met <- unique(plotdat$Metabolite)

for(i in 1:length(met)) {
  # create plot
  plot <- ggplot(data = filter(plotdat, Metabolite == met[i]), aes(x = treatment, y = ab.change, fill = treatment)) +
    geom_boxplot() +
    geom_jitter(width = 0.08, size = 2) +
    scale_fill_manual(values = treatmentIntcols) +
    theme_classic() +
    labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment") +
    ggtitle(as.character(met[i]))
  
  # save
  ggsave(filename = paste0("./data/plots/indiv-boxplot-plasmad35-ANOVA-", as.character(met[i]), ".png"), plot = plot, dpi = 600, height = 4.35, width = 7.32, units = "in")
  
  # print progress
  cat(as.character(met[i]))
}

# for the manuscript: need the number of replicates shown in each group
reps <- plotdat %>% 
  group_by(treatmentID, Metabolite) %>% 
  summarize(rep = length(unique(id)))

## PLOTS: EX-ER vs SED-ER (exercise-induced changes under energy restriction)
met <- c("D-4'-Phosphopantothenate", "Phenyllactic acid", "Xanthine")
exdat <- d35 %>% 
  filter(Metabolite %in% met) %>% 
  filter(!treatmentID %in% "SED-AL") 

# plot
ggplot(data = exdat, aes(x = treatmentID, y = ab.change, fill = treatmentID)) +
  geom_boxplot() +
  facet_wrap(~Metabolite) +
  scale_fill_manual(values = treatmentIntcols) +
  labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment")

## save 
#ggsave(filename = "./data/plots/changeSEDAL-plasma-EXERvsSEDER.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

### ---- KEGG functional groups -----

# get data
keg <- read.table("./data/aqeuousKEGG_Assignments.txt", sep = "\t", header = TRUE) %>% 
  # Subclass_4 does not have any metabolites; remove
  select(-Subclass_4) %>% 
  # create one column for each unique group
  mutate(allgroup = paste(Class, Subclass_1, Subclass_2, Subclass_3, sep = "_")) %>% 
  # remove "S" from beginning of mouse name
  mutate(id = str_remove_all(id, "S")) %>% 
  # parse down for joining
  select(Metabolite, allgroup) %>% 
  distinct()

# perform similar hypothesis testing

## ---- tumor KEGG groups ----

# get the average of each SED-AL metabolite

sedal <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # get only SED_AL
  filter(treatmentID == "SED-AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area),
            seSedAl = se(area))

# get change compared to SED-AL
abdat <- dat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  left_join(sedal, by = "Metabolite") %>% 
  # calculate change
  mutate(ab.change = area - avgSedAl) %>% 
  # SED-AL is now 0
  mutate(ab.change = case_when(
    treatmentID %in% "SED-AL" ~ 0,
    treatmentID != "SED-AL" ~ ab.change
  ))%>% 
  select(id, treatmentID, Metabolite, ab.change) %>% 
  # get functional annotations
  left_join(keg, by = "Metabolite") %>% 
  drop_na()

# perform ANOVA between treatment groups
group <- unique(abdat$allgroup)
pvals <- data.frame()

for(i in 1:length(group)) {
  
  # define model
  mod <- aov(ab.change ~ treatmentID, data = filter(abdat, allgroup == group[i]))
  
  # do Tukey posthoc
  tukey <- TukeyHSD(mod)
  
  # get p vals
  ps <- data.frame(Group = group[i],
                   contrast = rownames(tukey[1]$treatmentID),
                   pval = tukey[1]$treatmentID[,4],
                   row.names = NULL)
  
  # concatenate
  pvals <- rbind(pvals, ps)
  
}

# get significant 
sigs <- pvals %>% 
  filter(pval < 0.05)

# get table of mean +/- SE for results section
df <- abdat %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  group_by(treatmentID, Metabolite) %>% 
  summarize(avgchange = round(mean(ab.change), 3),
            sechange = round(se(ab.change), 3)) %>% 
  pivot_wider(names_from = treatmentID, values_from = c(avgchange, sechange))

# write to table
#write.table(df, "./data/changefromSEDAL-tumor-ANOVA.txt", sep = "\t", row.names = FALSE)

## PLOTS: Weight effects under exercise (EX-ER vs EX-AL)
wt <- c("5-Thymidylic acid", "Acetylphosphate", "ADP", "D-Glucose", "Hydroxproline", "Indole-3-carboxylic acid",
        "L-Alanine", "L-Arginine", "Succinic acid")
wtdat <- abdat %>% 
  filter(Metabolite %in% wt) %>% 
  filter(!treatmentID %in% "SED-AL")

ggplot(data = wtdat, aes(x = treatmentID, y = ab.change, fill = treatmentID)) +
  geom_boxplot() +
  facet_wrap(~Metabolite) +
  scale_fill_manual(values = treatmentIntcols) +
  labs(x = "Treatment", y = "Change from SED-AL", fill = "Treatment")

# save 
#ggsave(filename = "./data/plots/changeSEDAL-tumor-weighteffects.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

