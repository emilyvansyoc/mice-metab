## SUPPLEMENTARY FIGURES
# Emily Bean, 12/2020

# ---- get data ----

require(tidyverse)

# define standard error function
se <- function(x) sqrt(var(x)/length(x))

dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  # hyphenate instead of _ for treatment ID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# set ggplot theme 
theme_set(theme_classic())

# get colors
source("./R/RColorBrewer.R")

## get tissue type column for PCA
datf <- dat %>% 
  mutate(tissue.type = case_when(
    str_detect(id, "plasma") ~ "Plasma",
    str_detect(id, "tumor") ~ "Tumor"
  )) %>% 
  # make horizontal
  pivot_wider(names_from = "Metabolite", values_from = "area")

# matrix with only metab data
datpca <- datf %>% select(-c(id, treatmentID, Exercise, Weight, tissue.type))

## ---- S1: PCA of tissue type  ----

# make the PCA
pca <- PCA(datpca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph

# plot the biplot WITH LABELS 
fviz_pca_biplot(pca,
                # how to color the points
                col.in = datf$tissue.type,
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
ggsave(filename = "./data/plots/manuscript-plots/FigureS1.png", plot = last_plot(), dpi = 300)

### ---- S2A: PCA of plasma by treatment ----

# get just plasma
plasmadf <- datf %>% 
  filter(tissue.type == "Plasma")
plasmapca <- plasmadf %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, id))

# make the PCA
pca <- PCA(plasmapca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
pa <- fviz_pca_ind(pca,
             # how to color the points
             col.ind = plasmadf$treatmentID,
             legend.title = "Treatment",
             # shape of the points
             geom.ind = "point",
             addEllipses = TRUE,
             # only show the variables that have top 5 "contribution" to the pca
             #select.var = list(contrib = 10),
             #col.var = "black",
             #label = "none",
             # don't let them overlap
             #repel = TRUE,
             # x and y axis labels with variance
             xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
             ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
             title = "Plasma") +
  scale_color_manual(values = treatmentIntcols)


### ---- S2B: PCA of tumor by treatment ----

# get just tumor
tumordf <- datf %>% 
  filter(tissue.type == "Tumor")
tumorpca <- tumordf %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, id))

# make the PCA
pca <- PCA(tumorpca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
pb <- fviz_pca_ind(pca,
             # how to color the points
             col.ind = tumordf$treatmentID,
             legend.title = "Treatment",
             # shape of the points
             geom.ind = "point",
             addEllipses = TRUE,
             # only show the variables that have top 5 "contribution" to the pca
             #select.var = list(contrib = 10),
             #col.var = "black",
             #label = "none",
             # don't let them overlap
             #repel = TRUE,
             # x and y axis labels with variance
             xlab = paste0("PCA1 (", round(pca$eig[1, 2], 1), "%)"),
             ylab = paste0("PCA2 (", round(pca$eig[2, 2], 1), "%)"),
             title = "Tumor") +
  scale_color_manual(values = treatmentIntcols)

# arrange both and add letters
g1 <- ggarrange(pa, pb,
                ncol = 2, nrow = 1, common.legend = TRUE, legend = "right",
                labels = c("A", "B"))

# save
ggsave(filename = "./data/plots/manuscript-plots/FigureS2.png", plot = g1, dpi = 300)

### ---- S3: PCA of plasma by time ----

# get plasma dataframe and make Time column
pdf <- datf %>% 
  filter(tissue.type == "Plasma") %>% 
  mutate(Time = case_when(
    str_detect(id, "D7") ~ "Day 7",
    str_detect(id, "D21") ~ "Day 21",
    str_detect(id, "D35") ~ "Day 35"
  )) %>% 
  mutate(Time = factor(Time, ordered = TRUE, labels = c("Day 7", "Day 21", "Day 35")))
plasmapca <- pdf %>% 
  select(-c(treatmentID, Exercise, Weight, tissue.type, id, Time))

# make the PCA
pca <- PCA(plasmapca, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
fviz_pca_ind(pca,
                   # how to color the points
                   col.ind = pdf$Time,
                   legend.title = "Time",
                   # shape of the points
                   geom.ind = "point",
                   addEllipses = TRUE,
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
  scale_color_manual(values = timecols)
# save
ggsave(filename = "./data/plots/manuscript-plots/FigureS3.png", plot = last_plot(), dpi = 300)

