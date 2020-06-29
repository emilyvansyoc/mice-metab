
## PCA for Manuscript plots 

## ---- getData ----

# get data
source("./R/metabolomicsPCA.Rmd")

# set ggplot theme
theme_set(theme_minimal())

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

