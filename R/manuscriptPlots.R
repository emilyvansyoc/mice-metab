
## ---- PREP ----

# get colors
require(tidyverse)
require(ggpubr)
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

## ----- 10/2020 work ----
## FIG 1 is Sherry's experimental design

## ---- FIG. 2 Plasma corr to tumor volume ----

## Version 1: linear regression
# FIRST: run regression at ./R/aqueousTumorVol.R

# plot positive and negative slopes separately
pos <- sigs %>% filter(sign == "pos")
neg <- sigs %>% filter(sign == "neg")

plotdat <- df %>% 
  semi_join(neg, by = "Metabolite")

#ggplot(data = plotdat, aes(x = area, y = cm3)) +
 # geom_point(aes(shape = Time)) +
  #stat_smooth(se = FALSE, method = "lm", color = "black") +
  #facet_wrap(~Metabolite, scales = "free") +
  #labs(x = "Relative concentration", y = "Tumor volume cm3")

# create in ggpubr
p <- ggscatter(plotdat, x = "area", y = "cm3",
          add = "reg.line",
          conf.int = FALSE,
          shape = "Time",
          repel = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          add.params = list(color = "black",
                            fill = "black")) +
  facet_wrap(~ Metabolite, scales = "free") +
  #stat_cor(method = "pearson") +
  labs(x = "Relative concentration", y = expression(paste("Tumor volume ", cm^3))) 
# move the legend to the right
ggpar(p, legend = "right") 
# export the high quality plot _ FOR NOW JUST SCREENSHOT
#ggexport(g, filename = "./data/plots/manuscript-plots/Fig2-pos_slope.png", res = 300)




## ---- FIG 3: Plasma over time ----

# get plasma data from R/aqueousANOVAs.R

# get only PA+ER
paer <- sigs %>% filter(treatmentID == "PA+ER")
# filter data
plotdat <- abdat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(treatmentID == "PA+ER") %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 21", "Day 35")))

# Version 1: all 10 metabs on x axis, grouped by Time (may be too busy)
#ggplot(data = plotdat, aes(x = Metabolite, y = ab.change, fill = Time)) +
 # geom_boxplot() +
#  scale_fill_manual(values = timeGreys) +
 # labs(x = "Metabolite", y = "Change from SED+AL", fill = "Time") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Version 2: facet wrap by metabolite
#ggplot(data = plotdat, aes(x = Time, y = ab.change, fill = Time)) +
 # geom_boxplot() +
  #scale_fill_manual(values = timeGreys) +
  #labs(x = "Time", y = "Change from SED+AL", fill = "Time") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Metabolite, scales = "free")

## Version 2 in ggpubr
ggboxplot(data = plotdat, x = "Time", y = "ab.change", fill = "Time",
          facet.by = "Metabolite", #scales = "free",
          add = list("jitter", mean_sd),
          error.plot = "errorbar", 
          xlab = "Time", ylab = "\u0394 SED+AL") +
  scale_fill_manual(values = timeGreys) +
  theme(strip.background = element_rect(
     fill="white", linetype=0
  ))



## ---- FIG 4: Plasma at Day 35 ----


## Version 1: Y axis is relative concentration (use "plasma" dataframe)
## (run ANOVA from R/aqueousANOVAs.R before using dataframes)

# get only PA+ER vs SED+AL
paer <- sigs %>% filter(contrast == "SED+AL-PA+ER")

# filter data
plotdat <- d35 %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(!treatmentID == "SED+AL") %>% 
  mutate(treatmentID = factor(treatmentID, ordered = TRUE, levels = c("SED+AL", "PA+AL", "SED+ER", "PA+ER")))

# there's only 6 metabolites, so plot all in one 

#ggplot(data = plotdat, aes(x = treatmentID, y = area, fill = treatmentID)) +
 # geom_boxplot() +
#  scale_fill_manual(values = treatmentGreys) +
 # facet_wrap(~Metabolite, scales = "free") +
  #labs(x = "Treatment", y = "Relative concentration", fill = "Treatment") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Version 2: y axis is change SED+AL
# filter data
#plotdat <- d35 %>% 
 # semi_join(sigs, by = "Metabolite") %>% 
#  filter(!treatmentID == "SED+AL") %>% 
 # mutate(treatmentID = factor(treatmentID, ordered = TRUE, levels = c("PA+AL", "SED+ER", "PA+ER")))

# Version 2 in ggpubr
p <- ggboxplot(data = plotdat, x = "treatmentID", y = "ab.change", fill = "treatmentID",
          add = "jitter",
          xlab = "Treatment", ylab = "\u0394 SED+AL",
          title = "Xanthosine") +
  scale_fill_manual(values = treatmentGreys)
ggpar(p, legend = "none")

## ---- FIG. 5: Tumor corr to tumor volume ----

## Version 1: linear regression
#plotdat <- df %>% 
 # semi_join(sigs, by = "Metabolite")

#ggplot(data = plotdat, aes(x = area, y = cm3)) +
 # geom_point() +
  #stat_smooth(se = FALSE, method = "lm") +
  #facet_wrap(~Metabolite, scales = "free") +
  #labs(x = "Relative concentration", y = "Tumor volume cm3")

## Plot in ggpubr
# plot positive and negative slopes separately
pos <- sigs %>% filter(sign == "pos")
neg <- sigs %>% filter(sign == "neg")

plotdat <- df %>% 
  semi_join(neg, by = "Metabolite")

ggscatter(plotdat, x = "area", y = "cm3",
               add = "reg.line",
               conf.int = FALSE,
               repel = TRUE,
               cor.coef = TRUE,
               cor.method = "pearson",
               add.params = list(color = "black",
                                 fill = "black")) +
  facet_wrap(~ Metabolite, scales = "free") +
  #stat_cor(method = "pearson") +
  labs(x = "Relative concentration", y = expression(paste("Tumor volume ", cm^3))) +
  # make facet labels look better
  theme(strip.background = element_rect(
    fill="white", linetype=0
  ))



## ---- FIG. 6: Tumor differences between treatments (PA+ER and SED+AL) ----

# show all treatments with only the PA+ER vs SED+AL results
# (Similar to Fig 4)

# get only PA+ER vs SED+AL
paer <- sigs %>% filter(contrast == "SED+AL-PA+ER")

# filter data
plotdat <- abdat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(!treatmentID == "SED+AL") %>% 
  mutate(treatmentID = factor(treatmentID, ordered = TRUE, levels = c("SED+AL", "PA+AL", "SED+ER", "PA+ER")))

# plot in ggpubr
p <- ggboxplot(data = plotdat, x = "treatmentID", y = "ab.change", fill = "treatmentID",
               add = "jitter",
               xlab = "Treatment", ylab = "\u0394 SED+AL",
          facet.by = "Metabolite", scales = "free") +
  scale_fill_manual(values = treatmentGreys) +
  # make facet labels look better
  theme(strip.background = element_rect(
    fill="white", linetype=0
  ))
ggpar(p, legend = "none")

## ---- OLDER PLOTS ----
### ---- PCA ----
## PCA for Manuscript plots 

# get data
#source("./R/metabolomicsPCA.Rmd")


# get KEGG-cleaned metabolites
#source("./R/metabolomicsKEGG.Rmd") # the df we'll use is 'both'

# wrangle data
df <- both %>% 
  select(Metabolite = Match,
         id, treatmentID, Exercise, Weight, area)

# structure for PCA (needs to be horizontal)
dfh <- df %>% 
  pivot_wider(names_from = Metabolite, values_from = area) %>% 
  # get column for tissue type
  mutate(tissue.type = sapply(str_split(id, "_"), `[`, 2)) %>% 
  mutate(tissue.type1 = case_when(
    tissue.type %in% "tumor" ~ "Tumor",
    tissue.type %in% "plasmaD7" ~ "Plasma",
    tissue.type %in% "plasmaD21" ~ "Plasma",
    tissue.type %in% "plasmaD35" ~ "Plasma"
  )) %>% 
  column_to_rownames(var = "id")

# remove extra columns for PCA
dfpca <- dfh %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, tissue.type1))

## ---- PCA by Tissue Type ----

# make the PCA
pca <- PCA(dfpca, 
            scale.unit = TRUE, # scale everything to equal variance
            graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
fviz_pca_biplot(pca,
                # how to color the points
                col.in = dfh$tissue.type1,
                legend.title = "Tissue Type",
                # shape of the points
                geom.ind = "point",
                # no transparency
                alpha.ind = 1,
                # add ellipses at 95% confidence level
                addEllipses = TRUE,
                ellipse.level = 0.95,
                # only show the variables that have top 5 "contribution" to the pca
                select.var = list(contrib = 10),
                col.var = "black",
                #label = "none",
                # don't let them overlap
                repel = TRUE,
                # x and y axis labels with variance
                xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
                ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
                title = NULL)+
  scale_color_manual(values = tissuecols)

# save
ggsave(filename = "./data/plots/pca-plasma-tumor-labels.tiff", plot = last_plot(), dpi = "print")

## PCA with NO LABELS 

fviz_pca_biplot(pca,
                # how to color the points
                col.in = dfh$tissue.type1,
                legend.title = "Tissue Type",
                # shape of the points
                geom.ind = "point",
                # no transparency
                alpha.ind = 1,
                # add ellipses at 95% confidence level
                addEllipses = TRUE,
                ellipse.level = 0.95,
                # only show the variables that have top 5 "contribution" to the pca
                select.var = list(contrib = 10),
                col.var = "black",
                label = "none",
                # don't let them overlap
                repel = TRUE,
                # x and y axis labels with variance
                xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
                ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
                title = NULL)+
  scale_color_manual(values = tissuecols)
# save
ggsave(filename = "./data/plots/pca-plasma-tumor-nolabels.tiff", plot = last_plot(), dpi = "print")

## ---- PCA 4 treatment groups: Tumor ----

# get just tumor
tumordf <- dfh %>% 
  filter(tissue.type1 == "Tumor")
tumorpca <- tumordf %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, tissue.type1))

# make the PCA
pca <- PCA(tumorpca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
fviz_pca_ind(pca,
                # how to color the points
                col.ind = tumordf$treatmentID,
                legend.title = "Treatment",
                # shape of the points
                geom.ind = "point",
                #addEllipses = TRUE,
                # only show the variables that have top 5 "contribution" to the pca
                #select.var = list(contrib = 10),
                #col.var = "black",
                #label = "none",
                # don't let them overlap
                #repel = TRUE,
                # x and y axis labels with variance
                xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
                ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
                title = NULL) +
  scale_color_manual(values = treatmentIntcols)

# save
ggsave(filename = "./data/plots/pca-tumor-treatments.tiff", plot = last_plot(), dpi = "print")

## ---- PCA 4 treatment groups: Plasma ----

# get just tumor
plasmadf <- dfh %>% 
  filter(tissue.type1 == "Plasma")
plasmapca <- plasmadf %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, tissue.type1))

# make the PCA
pca <- PCA(plasmapca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
fviz_pca_ind(pca,
             # how to color the points
             col.ind = plasmadf$treatmentID,
             legend.title = "Treatment",
             # shape of the points
             geom.ind = "point",
             #addEllipses = TRUE,
             # only show the variables that have top 5 "contribution" to the pca
             #select.var = list(contrib = 10),
             #col.var = "black",
             #label = "none",
             # don't let them overlap
             #repel = TRUE,
             # x and y axis labels with variance
             xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
             ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
             title = NULL) +
  scale_color_manual(values = treatmentIntcols)

# save
ggsave(filename = "./data/plots/pca-plasma-treatments.tiff", plot = last_plot(), dpi = "print")

## ---- KEGG classes ----

#source("./R/metabolomicsKEGG.Rmd")

## replace _ with - for plots
assignv <- assignv %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

theme_set(theme_minimal())

## tumor - barplot of metabolic classes
ggplot(data = filter(assignv, tissue.type == "tumor"), 
       aes(x = Class, y = area, fill = treatmentID)) +
  geom_col() +
  scale_fill_manual(values = treatmentIntcols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Metabolic Class", y = "Area Under the Curve", fill = "Treatment")

#ggsave(filename = "./data/plots/barplot-tumor-treatment-class.tiff", plot = last_plot(), dpi = "print")

## plasma - barplot of metabolic classes
ggplot(data = filter(assignv, tissue.type == "plasma"), 
       aes(x = Class, y = area, fill = treatmentID)) +
  geom_col() +
  scale_fill_manual(values = treatmentIntcols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Metabolic Class", y = "Area Under the Curve", fill = "Treatment")

#ggsave(filename = "./data/plots/barplot-plasma-treatment-class.tiff", plot = last_plot(), dpi = "print")

### ----- ANOVA barplots - OLD ----

#source("./R/aqueousANOVAs.R")

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

## TUMOR

# get significant interaction terms
df <- as.data.frame(tumormod[2]) %>% 
  filter(comparison == "Interaction")
df # Hydroxyproline.Aminolevulinate, Quinolinate, Acetyl.aspartate, Glucose

# get these metabolites
mets <- c("Hydroxyproline.Aminolevulinate", "Quinolinate", "Acetyl.aspartate", "Glucose")
subdf <- aq %>% 
  filter(tissue.type == "tumor") %>% 
  filter(metabolite %in% mets) %>% 
  # get the correct metabolite names 
  inner_join(names, by = c("metabolite" = "oldnames")) %>% 
  mutate(Metabolite = case_when(
    is.na(KEGGMatch) ~ cleanednames,
    !is.na(KEGGMatch) ~ KEGGMatch
  )) %>% 
  # calculate mean & se
  group_by(Metabolite, treatmentID) %>% 
  summarize(meanarea = mean(area),
            searea = se(area)) %>% 
  ungroup() %>% 
  # change _ to - in treatmentID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# make plot
ggplot(data = subdf, aes(x = treatmentID, y = meanarea, fill = Metabolite)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = meanarea - searea, ymax = meanarea + searea), 
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_manual(values = accentcols) %>% 
  labs(x = "Treatment", y = "Area Under the Curve") +
  facet_wrap(~Metabolite)
#ggsave(filename = "./data/plots/barplot-tumor-ANOVA-interactions-facetwrap.tiff", plot = last_plot(), dpi = "print")

### PLASMA

pdf <- as.data.frame(plasmamod[2]) %>% 
  # get only interactions
  filter(comparison == "Ex-Wt-Time")

# subset data to get significant interactions
psub <- aq %>% 
  filter(tissue.type == "plasma") %>% 
  left_join(names, by = c("metabolite" = "oldnames")) %>% 
  select(-metabolite.y) %>% 
  mutate(Metab.name = case_when(
    !is.na(Match) ~ Match,
    is.na(Match) ~ metabolite
  )) %>% 
  # get only some metabs
  semi_join(pdf, by = "metabolite") %>% 
  
  # change column names to only keep correct metabolite names
  rename(Metabolite = Metab.name) %>% 
  # calculate mean & se
  group_by(Metabolite, treatmentID, Label) %>% 
  summarize(meanarea = mean(area),
            searea = se(area)) %>% 
  ungroup() %>% 
  # change _ to - in treatmentID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"),
         Time = case_when(
           Label %in% "plasmaD7" ~ "Day 7",
           Label %in% "plasmaD21" ~ "Day 21",
           Label %in% "plasmaD35" ~ "Day 35"
         )) %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 21", "Day 35")))

ggplot(data = psub, aes(x = treatmentID, y = meanarea, fill = Metabolite)) +
  geom_col(position = position_dodge(0.9)) +
  facet_wrap(~Time) +
  labs(x = "Treatment", y = "Scaled Concentration")
#ggsave(filename = "./data/plots/barplot-plasma-ANOVA-interactions.tiff", plot = last_plot(), dpi = "print")


## ---- Tumor tissue regressions ----

# source("./R/metabTumorWeightRegression.Rmd")

# run just Tumor chunk

# get correct metabolite names

tumdf <- sigdf %>% 
  left_join(names, by = c("metabolite" = "oldnames")) %>% 
  mutate(Metabolite = case_when(
    is.na(KEGGMatch) ~ cleanednames,
    !is.na(KEGGMatch) ~ KEGGMatch
  ))

ggplot(data = tumdf, aes(x = area, y = cm3)) +
  geom_point() +
  # add best fit line without 95% confidence interval
  geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, fullrange = TRUE) +
  facet_wrap(~Metabolite, scales = "free") +
  labs(x = "Scaled area", y = "Tumor volume (cm3)")
#ggsave(filename = "./data/plots/tumor-reg-tumortissue.tiff", plot = last_plot(), dpi = "print")

## run just PLASMA chunk

# now there is a Time consideration 

plasdf <- sigdf %>% 
  left_join(names, by = c("metabolite" = "oldnames")) %>% 
  # not all metabolites have a KEGG match
  mutate(Metabolite = case_when(
    is.na(KEGGMatch) ~ cleanednames,
    !is.na(KEGGMatch) ~ KEGGMatch
  ),
  Time = factor(Time, ordered = TRUE, levels = c("D7", "D21", "D35")))

# plot with Time as facet grid
ggplot(data = plasdf, aes(x = area, y = cm3)) +
  geom_point() +
  # add best fit line without 95% confidence interval
  geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, fullrange = TRUE) +
  facet_wrap(Time ~ Metabolite, scales = "free") +
  labs(x = "Scaled area", y = "Tumor volume (cm3)")
#ggsave(filename = "./data/plots/tumor-reg-plasma.tiff", plot = last_plot(), dpi = "print")


  
