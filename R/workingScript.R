
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

BiocManager::install("ComplexHeatmap")
require(ComplexHeatmap)

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
