
## PCA for metabolites
# FIRST - use transformed data... will this work? 

## ALL - TISSUE TYPE

# get only variables, in horizontal form
aqh <- aq %>% 
  spread(key = metabolite, value = area)
# perform PCA
pca <- prcomp(aqh[, 8:ncol(aqh)], center = TRUE, scale. = TRUE)
summary(pca) # PCA 1 and 2 explain most of the variance


library(devtools)
install_github("vqv/ggbiplot")
require(ggbiplot)

ggbiplot(pca, var.axes = FALSE, groups = aqh$tissue.type, ellipse = TRUE) +
  ggtitle("Tissue type")


# perform a PCA for just tumor tissue

tumorh <- aqh %>% filter(tissue.type == "tumor")

pcatumor <- prcomp(tumorh[, 8:ncol(tumorh)], center = TRUE, scale. = TRUE)
summary(pcatumor)

# biplot
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$treatmentID, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Weight, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Exercise, ellipse = TRUE)

## do these look better on untransformed data?
# perform PCA
pca <- prcomp(aqproch[, 6:ncol(aqproch)], center = TRUE, scale. = TRUE)
summary(pca)
ggbiplot(pca, var.axes = FALSE, groups = aqproch$tissue.type, ellipse = TRUE)

# TUMOR
pcatumor <- prcomp(tumorh[, 6:ncol(tumorh)], center = TRUE, scale. = TRUE)
summary(pcatumor)

# biplot
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$treatmentID, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Weight, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Exercise, ellipse = TRUE)

## PLASMA
pcaplasma <- prcomp(plasmaproch[, 7:ncol(plasmaproch)], center = TRUE, scale. = TRUE)
# error - there are columns with all zeros/constants
# did Metaboanalyst replace zeros with 358.5??



summary(pcatumor)

# biplot
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$treatmentID, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Weight, ellipse = TRUE)
ggbiplot(pcatumor, var.axes = FALSE, groups = tumorh$Exercise, ellipse = TRUE)
