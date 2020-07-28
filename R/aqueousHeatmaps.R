## Aqueous: heatmaps

## ---- get data ----

require(ComplexHeatmap)
require(tidyverse)

source("./R/RColorBrewer.R")

dat <- read.table("./data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# ---- tumor pre-process -----

# get the average of each SED-AL metabolite

sedal <- dat %>% 
  # get only tumor
  filter(str_detect(id, "_tumor")) %>% 
  # get only SED_AL
  filter(treatmentID == "SED-AL") %>% 
  group_by(Metabolite) %>% 
  # average
  summarize(avgSedAl = mean(area))

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

# heatmaps need to be in matrix format
df <- abdat %>% 
  filter(!treatmentID == "SED-AL") %>% 
  mutate(newid = str_extract_all(id, "M(\\d){1,2}")) %>% 
  select(newid, treatmentID, ab.change, Metabolite) %>% 
  pivot_wider(names_from = Metabolite, values_from  = ab.change) %>% 
  column_to_rownames(var = "newid") 
dfmat <- df %>% select(-treatmentID)
dfmat <- as.matrix(dfmat)

# plot heatmap of all metabolites
graphics.off()
#pdf("./data/plots/heatmap-tumor-allmetabolites.pdf", width = 11, height = 11)
Heatmap(t(dfmat), row_dend_reorder = TRUE, use_raster = FALSE, top_annotation = HeatmapAnnotation(Treatment = df$treatmentID, 
                                                       col = list(Treatment = c("EX-AL" = "#FDBF6F", "EX-ER" = "#FF7F00",  "SED-ER" = "#6A3D9A" ))))
dev.off()

## Heatmap of just significant metabolites (from ./R/aqueousANOVAs.R)
mets <- c("5-Thymidylic acid", "Acetylphosphate", "ADP", "D-Glucose", "Hydroxyproline", "Indole-3-carboxylic acid", "L-Alanine",
          "L-Arginine", "L-Cystathionine", "L-Glutamine", "L-Histidine", "L-Lysine", "N-Acetyl-L-aspartic acid",
          "Quinolinic acid", "Succinic acid", "Taurine", "Uridine triphosphate")

df1 <- abdat %>% 
  filter(!treatmentID == "SED-AL") %>% 
  filter(Metabolite %in% mets) %>% 
  mutate(newid = str_extract_all(id, "M(\\d){1,2}")) %>% 
  select(newid, treatmentID, ab.change, Metabolite) %>% 
  pivot_wider(names_from = Metabolite, values_from  = ab.change) %>% 
  column_to_rownames(var = "newid") 
dfmat1 <- df1 %>% select(-treatmentID)
dfmat1 <- as.matrix(dfmat1)

#pdf("./data/plots/heatmap-tumor-sigmetabolites.pdf", width = 11, height = 11)
Heatmap(t(dfmat1), row_dend_reorder = TRUE, use_raster = FALSE, top_annotation = HeatmapAnnotation(Treatment = df$treatmentID, 
                                                                                                  col = list(Treatment = c("EX-AL" = "#FDBF6F", "EX-ER" = "#FF7F00",  "SED-ER" = "#6A3D9A" ))))
dev.off()
