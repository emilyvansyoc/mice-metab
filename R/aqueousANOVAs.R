
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

## ---- PLASMA ----

# get only plasma and calculate Time column
plasma <- dat %>% 
  filter(str_detect(id, "_plasma")) %>% 
  # make Time column
  mutate(Time = case_when(
    str_detect(id, "_plasmaD7") ~ "Day 7",
    str_detect(id, "_plasmaD21") ~ "Day 21",
    str_detect(id, "_plasmaD35") ~ "Day 35"
  ))


# get only D7
d7 <- plasma %>% 
  filter(Time == "Day 7") %>% 
  group_by(Metabolite) %>% 
  summarize(avgD7 = mean(area),
            seD7 = se(area))

# calculate change from D7
tdat <- plasma %>% 
  left_join(d7, by = "Metabolite") %>% 
  mutate(ab.change = area - avgD7) %>% 
  # change D7 to zero
  mutate(ab.change = case_when(
    Time == "Day 7" ~ 0,
    Time != "Day 7" ~ ab.change
  )) %>% 
  select(id, treatmentID, Metabolite, Time, ab.change) %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 21", "Day 35")))

# perform ANOVA for each metabolite and treatment to see changes over time
metab <- unique(tdat$Metabolite)
trt <- unique(tdat$treatmentID)
pvals <- data.frame()

for(j in 1:length(metab)) {
  
  for(i in 1:length(trt)) {
    
    # define model
    mod <- aov(ab.change ~ Time, data = filter(tdat, Metabolite == metab[j] & treatmentID == trt[i]))
    
    # do post-hoc
    tukey <- TukeyHSD(mod)
    
    # get p vals
    ps <- data.frame(Metabolite = metab[j],
                     Treatment = trt[i],
                     contrast = rownames(tukey[1]$Time),
                     pval = round(tukey[1]$Time[,4], 3),
                     row.names = NULL)
    
    # concatenate
    pvals <- rbind(pvals, ps)
    
  }
}

# get significance
sigs <- pvals %>% filter(pval < 0.05)

sigEXER <- sigs %>% filter(Treatment == "EX-ER")
# what's a good way to visualize this?
sigdf <- tdat %>% 
  semi_join(sigEXER, by = "Metabolite") %>% 
  group_by(Metabolite, treatmentID, Time) %>% 
  summarize(avg = mean(ab.change),
            sd = se(ab.change))

ggplot(data = filter(sigdf, treatmentID == "EX-ER"), aes(x = Time, y = avg, group = Metabolite)) +
  geom_point() +
  geom_line()
