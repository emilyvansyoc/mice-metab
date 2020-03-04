# read all lipids data
lip <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/allLipidsCleaned.txt", header = TRUE, stringsAsFactors = FALSE) 

## Make all data horizontal for PCAs
liph <- lip %>% 
  select(-ion) %>% 
  spread(key = lipid, value = log.area)

# plasma only
plasmah <- liph %>% 
  filter(tissue.type == "plasma") %>% 
  mutate(Time = sapply(strsplit(SampleName, "_"), `[`, 1)) %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("D7", "D21", "D35")))

# tumor only
tumorh <- liph %>% 
  filter(tissue.type == "tumor")

# save Rdata to load faster
#save.image("./data/parsedLipids.RData")