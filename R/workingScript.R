# anovas and kegg visualizations

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

## test plot
test <- assignv %>% 
  # make a smaller dataset
  slice_sample(n = 10)

ggplot(data = test, aes(x = Exercise, y = area, fill = Metabolite)) +
  geom_col()

# get significant data
se <- function(x) sqrt(var(x)/length(x))

# get significant interaction terms
df <- as.data.frame(tumormod[2])
df # Hydroxyproline.Aminolevulinate, Quinolinate, Acetyl.aspartate, Glucose

# get these metabolites
mets <- c("Hydroxyproline.Aminolevulinate", "Quinolinate", "Acetyl.aspartate", "Glucose")
subdf <- aq %>% 
  filter(tissue.type == "tumor") %>% 
  filter(metabolite %in% mets) %>% 
  # get the correct metabolite names 
  inner_join(names, by = c("metabolite" = "oldnames")) %>% 
  # change column names to only keep correct metabolite names
  select(-metabolite, Metabolite = Match) %>% 
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
  labs(x = "Treatment", y = "Area Under the Curve")
#ggsave(filename = "./data/plots/barplot-tumor-ANOVA-interactions.tiff", plot = last_plot(), dpi = "print")

# test PCA on functional groups

test <- assignv %>% 
  filter(tissue.type == "tumor") %>% 
  group_by(Class, id, treatmentID) %>% 
  summarize(meanarea = mean(area)) %>% 
  pivot_wider(names_from = Class, values_from = meanarea) %>% 
  column_to_rownames(var = "id") %>% 
  select(-treatmentID)

testh <- assignv %>% 
  filter(tissue.type == "tumor") %>% 
  group_by(Class, id, treatmentID) %>% 
  summarize(meanarea = mean(area)) %>% 
  pivot_wider(names_from = Class, values_from = meanarea) %>% 
  column_to_rownames(var = "id") 

# make the PCA
pca <- PCA(test, 
           scale.unit = TRUE, # scale everything to equal variance
           graph = FALSE) # don't print a graph


# plot the biplot WITH LABELS 
fviz_pca_biplot(pca,
                # how to color the points
                col.in = testh$treatmentID,
                legend.title = "Tissue Type",
                # shape of the points
                geom.ind = "point",
                # no transparency
                alpha.ind = 1,
                # add ellipses at 95% confidence level
                #addEllipses = TRUE,
                #ellipse.level = 0.95,
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
  scale_color_manual(values = treatmentIntcols)

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


# get subset data with Classes to categorize
plasdf <- assignv %>% 
  select(Metabolite, Class, Subclass_1, treatmentID) %>% 
  distinct() %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-")) %>% 
  # get all metabolites, including the 4 without KEGG identifiers
  right_join(psub, by = c("Metabolite", "treatmentID")) %>% 
  # change NA to "No category" 
  mutate(Class = as.character(Class)) %>% 
  mutate(Class = case_when(
    is.na(Class) ~ "No Category",
    Class %in% "Vitamins and Cofactors" ~ "Vitamins",
    Class %in% "Phytochemical compounds" ~ "Phytochemicals",
    !is.na(Class) ~ Class
  ))
  


# Plot everything in a facet grid
ggplot(data = plasdf, aes(x = treatmentID, y = meanarea, fill = Metabolite))+
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanarea - searea, ymax = meanarea + searea), 
                position = position_dodge(width = 0.9), width = 0.2) +
  facet_grid(Class~Time, scales = 'free') + 
  scale_fill_manual(values = getPalette(length(unique(test$Metabolite)))) +
  labs(x = "Treatment", y = "Scaled Concentration")
#ggsave(filename = "./data/plots/barplot-facetgrid-plasma-ANOVA-interactions.tiff", plot = last_plot(), dpi = "print")