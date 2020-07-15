
## Summary statistics for interpreting ANOVA results

# source("~/R/metabolomicsKEGG.Rmd")

sum <- plasma %>% 
  group_by(treatmentID, metabolite, Time) %>% 
  summarize(meanarea = mean(area),
            searea = se(area)) %>% 
  ungroup() %>% 
  left_join(names, by = c("metabolite" = "oldnames")) %>% 
  mutate(Metabolite = case_when(
    is.na(KEGGMatch) ~ cleanednames,
    !is.na(KEGGMatch) ~ KEGGMatch
  )) %>% 
  select(-c(cleanednames, KEGGMatch, metabolite))

# for saving; spread by treatmentID
psp <- sum %>% 
  pivot_wider(names_from = c(treatmentID, Time), values_from = c(meanarea, searea))

#write.table(psp, file = "./data/sumstats-plasma-cleanednames.txt", sep = "\t", row.names = FALSE)

# Heatmaps

#BiocManager::install("ComplexHeatmap")
require(ComplexHeatmap)

mets <- c("Kynurenic acid", "Quinolinic acid", "L-Tryptophan")

dat <- read.table("./data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  filter(str_detect(id, "_tumor")) %>% 
  filter(Metabolite %in% mets) %>% 
  pivot_wider(names_from = Metabolite, values_from = area) %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))
subdat <- dat %>% 
  select(-c(treatmentID, Exercise, Weight)) %>% 
  column_to_rownames(var = "id")
subdat <- as.matrix(subdat)


Heatmap(t(subdatf), top_annotation = HeatmapAnnotation(Treatment = datf$treatmentID, 
                                                      col = list(Treatment = c("EX-AL" = "#FDBF6F", "EX-ER" = "#FF7F00", "SED-AL" = "#CAB2D6", "SED-ER" = "#6A3D9A" ))))

# get only treatment vs control
trt <- c("EX-ER", "SED-AL")
datf <- dat %>% 
  filter(treatmentID %in% trt)
subdatf <- datf %>% 
  select(-c(treatmentID, Exercise, Weight)) %>% 
  column_to_rownames(var = "id")

# tumor
tumh <- tumor %>% 
  left_join(names, by = c("metabolite" = "oldnames")) %>% 
  mutate(Metabolite = case_when(
    is.na(KEGGMatch) ~ cleanednames,
    !is.na(KEGGMatch) ~ KEGGMatch
  )) %>% 
  select(-c(cleanednames, KEGGMatch, metabolite)) %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-")) %>% 
  # make into matrix for heatmap 
  mutate(rowid = seq(1:nrow(tumh))) %>% 
  select(id, Metabolite, area) %>% 
  pivot_wider(names_from = Metabolite, values_from = area)

## metaboanalyst
trt <- c("EX_ER", "SED_AL", sep = "|")
dat <- read.table("./data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE)
df <- dat %>% 
  select(id, treatmentID, area, Metabolite) %>% 
  filter(str_detect(id, "_tumor")) %>% 
  pivot_wider(names_from = Metabolite, values_from = area) %>% 
  filter(treatmentID %in% trt) %>% 
  mutate(treatment = case_when(
    treatmentID %in% "EX_ER" ~ "Treatment",
    treatmentID %in% "SED_AL" ~ "Control"
  )) %>% 
  select(-treatmentID) %>% 
  select(id, treatment, `Pyruvic acid`:ncol(df))

write.table(df, "./data/tumorMetaboAnalyst.txt", sep = "\t", row.names = FALSE)


# get un-normalized data, only tumor
dat <- read.table("./data/aqueousCleanedNotNormalized.txt", sep = "", header = TRUE)
df <- dat %>% 
  filter(Label == "tumor") %>% 
  select(-c(id, Exercise, Weight)) %>% 
  filter(treatmentID %in% trt) %>% 
  mutate(treatment = case_when(
    treatmentID %in% "EX_ER" ~ "Treatment",
    treatmentID %in% "SED_AL" ~ "Control"
  )) %>% 
  select(-treatmentID) %>% 
  pivot_wider(names_from = metabolite, values_from = area) %>% 
  select(-Label)
write.table(df, "./data/tumorNotNormalized.txt", sep = "\t", row.names = FALSE)


## what is going on with M1
datv <- datf %>% 
  pivot_longer(cols = 5:ncol(datf), names_to = "Metabolite", values_to = "area") %>% 
  mutate(id = str_replace_all(id, "_tumor", ""))

ggplot(data = datv, aes(x = id, y = area, fill = treatmentID)) +
  geom_col() +
  facet_wrap(~Metabolite) +
  theme(axis.text.x = element_text(angle = 45,  hjust=1))


## ---- PERCENT CHANGE FROM SED_AL ----
dat <- read.table("./data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE)
pdat <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  select(-c(Exercise, Weight)) %>% 
  group_by(Metabolite, treatmentID) %>% 
  summarize(mean = mean(area)) %>% 
  pivot_wider(names_from = treatmentID, values_from = mean) %>% 
  pivot_longer(cols = c(EX_AL, EX_ER, SED_ER), names_to = "Treatment", values_to = "avg") %>% 
  mutate(percent.change = ((avg - SED_AL) / avg) * 100) %>% 
  select(-c(SED_AL, avg)) %>% 
  pivot_wider(names_from = Metabolite, values_from = percent.change)

ggplot(data = pdat, aes(x =  Treatment, y = `Quinolinic acid`, fill = Treatment)) +
  geom_col()

# try by average SED_AL and then calculating for each mouse

# get averages for SED_AL
se <- function(x) sqrt(var(x)/length(x))
sedal <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # get only SED_AL
  filter(treatmentID == "SED_AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area),
            seSedAl = se(area)) %>% 
  # make dummy column for Metabolite for plot below
  mutate(metabId = seq(1:nrow(sedal))) 
 
# diagnostic: what is the variance in SEDAL?
ggplot(data = sedal, aes(x = metabId, y = avgSedAl)) +
  geom_col() +
  geom_errorbar(aes(ymin = avgSedAl - seSedAl, ymax = avgSedAl + seSedAl)) +
  theme_bw() +
  labs(x = "Metabolite Index", y = "Mean +/- Standard Error", 
       title = "SED-AL Metabolites have wide variance")

# get percent change compared to SED-AL
newdat <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # remove SED-AL
  #filter(!treatmentID == "SED_AL") %>% 
  # combine with SED-AL averages
  left_join(sedal, by = "Metabolite") %>% 
  # calculate percent change
  mutate(percent.change = ((area - avgSedAl) / area) * 100) %>% 
  mutate(percent.change = case_when(
    treatmentID %in% "SED_AL" ~ 0,
    treatmentID != "SED_AL" ~ percent.change
  )) %>% 
  select(id, treatmentID, Metabolite, percent.change) 

# what about just absolute change without percent
abdat <- dat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  left_join(sedal, by = "Metabolite") %>% 
  mutate(ab.change = area - avgSedAl) %>% 
  mutate(ab.change = case_when(
    treatmentID %in% "SED_AL" ~ 0,
    treatmentID != "SED_AL" ~ ab.change
  ))%>% 
  select(id, treatmentID, Metabolite, ab.change) 
  
# calculate fold change
fdat <- dat %>% 
  filter(str_detect(id, "_tumor")) %>% 
  left_join(sedal, by = "Metabolite") %>% 
  mutate(f.change = area / avgSedAl) %>% 
  mutate(f.change = case_when(
    treatmentID %in% "SED_AL" ~ 0,
    treatmentID != "SED_AL" ~ f.change
  ))%>% 
  select(id, treatmentID, Metabolite, f.change) 

ggplot(data = newdat, aes(x = treatmentID, y = `Quinolinic acid`, fill = treatmentID)) +
  geom_col()

## ANOVA with new percent change structure

metabs <- unique(fdat$Metabolite)
pvals <- data.frame()

for(i in 1:length(metabs)) {
  
  # define model
  mod <- aov(f.change ~ treatmentID, data = filter(fdat, Metabolite == metabs[i]))
  
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

# significant 
sigs <- pvals %>% 
  filter(pval < 0.05)

# plot one to look at
ggplot(data = filter(fdat, Metabolite == "Quinolinic acid"), aes(x = treatmentID, y = f.change)) +
  geom_boxplot()

# make paneled plot of all significant metabolites
sigdf <- fdat %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  filter(!treatmentID == "SED_AL")

# plot all
ggplot(data = sigdf, aes(x = treatmentID, y = f.change)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~Metabolite)
