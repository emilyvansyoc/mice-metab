## PCA in ggfortify

require(ggfortify)


#  PCA

# Prepare the data
df <- iris[, -5]

aqh <- aq %>% pivot_wider(names_from = metabolite, values_from = area) 
dat <- aqh %>% select(-c(id, Label, mouseID, treatmentID, Exercise, Weight, tissue.type))
# Principal component analysis
pca <- prcomp(dat, scale. = TRUE)
# Plot
autoplot(pca, loadings = FALSE, loadings.label = FALSE,
         data = aqh, colour = 'treatmentID')
