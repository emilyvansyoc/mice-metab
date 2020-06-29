## Aqueous samples multivariate analysis

## ---- readData ----

# this is processed and normalized data from MetaboAnalyst
require(dplyr)
require(tidyr)
require(ggplot2)
aq <- read.csv("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/data_normalizedAQ.csv", header = TRUE, stringsAsFactors = TRUE) %>% 
  mutate(id = X,
         X = NULL)
# read in processed but NON-normalized data from MetaboAnalyst for community analysis
aqproc <- read.csv("./data/data_processedAQ.csv", header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(id = X)

# add sample keys
key <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/sampleKey.txt", header = TRUE) %>% 
  # add runID with a dash
  mutate(tumor = paste(runID, tumor, sep = "-"),
         plasmaD7 = paste(runID, plasmaD7, sep = "-"),
         plasmaD21 = paste(runID, plasmaD21, sep = "-"),
         plasmaD35 = paste(runID, plasmaD35, sep = "-")) %>% 
  # make vertical
  gather(key = "Label", value = "id", tumor, plasmaD7, plasmaD21, plasmaD35) %>% 
  mutate(runID = NULL)


# merge to join IDs

all <- merge(key, aqproc, by = c("id", "Label")) %>% 
  mutate(treatmentID = as.character(treatmentID))

for(i in 1:nrow(all)) {
  all$Exercise[i] <- strsplit(all$treatmentID[i], "_")[[1]][1]
  all$Weight[i] <- strsplit(all$treatmentID[i], "_")[[1]][2]
}

# reorder to look better
aqdata <- all[, c(1:4, 141:142, 5:140)] 

## make vertical copy to work with
aqdataVert <- aqdata %>% 
  gather(key = "metabolite", value = "area", 7:ncol(aqdata)) %>% 
  mutate(id = factor(id),
         Label = factor(Label),
         mouseID = factor(mouseID),
         treatmentID = factor(treatmentID),
         Exercise = factor(Exercise),
         Weight = factor(Weight),
         metabolite = factor(metabolite))

# double check - are we meeting assumptions of normality?
hist(aqdataVert$area)
qqnorm(aqdataVert$area)
qqline(aqdataVert$area)

# write for RMD
#write.table(aqdataVert, file = "./data/aqueousCleanedNotNormalized.txt", row.names = FALSE)

