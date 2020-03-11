
# Building a heatmap 

df <- tukdf %>% 
  spread(key = metabolite, value = area) %>% 
  select(Acetoacetate:ncol(df))
df <- as.matrix(df)

heatmap(t(df))
require(tibble)
sum <- tukdf %>% 
  group_by(Exercise, Weight, treatmentID, metabolite) %>% 
  summarize(meanarea = mean(area)) %>% 
  ungroup() %>% 
  spread(key = metabolite, value = meanarea) %>% 
  column_to_rownames(var = "treatmentID") %>% 
  select(Acetoacetate:ncol(sum))



sum <- as.matrix(sum)

heatmap(t(sum))
