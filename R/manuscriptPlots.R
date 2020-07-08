
## ---- PREP ----

# get colors
source("./R/RColorBrewer.R")

# set ggplot theme
theme_set(theme_minimal())

# function to calculate standard error
se <- function(x) sqrt(var(x)/length(x))

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

## ---- Tissue Type ----

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

## ---- 4 treatment groups: Tumor ----

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

## ---- 4 treatment groups: Plasma ----

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

### ----- ANOVA barplots ----

#source("./R/metabolomicsANOVAs.Rmd")

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
