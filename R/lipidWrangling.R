## LIPIDS: Data Wrangling
require(dplyr)
require(tidyr)
require(stringr)

sampkey <- read.csv("lipid_SampleKeys1.csv", stringsAsFactors = FALSE)[-1,] # remove Blank row

mousekey <- read.table("mouseKeys.txt", header = TRUE, stringsAsFactors = FALSE)


## wrangle the heck out of this

for(i in 1:nrow(mousekey)) {
  mousekey[i, "id"] <- strsplit(mousekey$mouseID[i], "-")[[1]][3]
  mousekey[i, "Exercise"] <- strsplit(mousekey$treatmentID[i], "_")[[1]][1]
  mousekey[i, "Weight"] <- strsplit(mousekey$treatmentID[i], "_")[[1]][2]
}

for(i in 1:nrow(sampkey)) {
  
  # not all mouseIDs have length 2
  if(grepl("_", sampkey$SampleName[i]) == TRUE) {
    
    # these are tumor samples; kill 2 birds w/1 stone
    sampkey[i, "id"] <- strsplit(sampkey$SampleName[i], "_")[[1]][2]
    sampkey[i, "tissue.type"] <- "tumor"
    
  } else {
    sampkey[i, "id"] <- sampkey$SampleName[i]
    sampkey[i, "tissue.type"] <- "plasma"
  }
  
  # separate positive and negative
  if(grepl("POS", sampkey$SampleID[i]) == TRUE) {
    
    sampkey[i, "ion"] <- "POS"
  } else {
    sampkey[i, "ion"] <- "NEG"
  }
  
  # for all samples: remove .wiff to identify w/metaboanalyst CSV
  sampkey[i, "shortID"] <- strsplit(sampkey$SampleID[i], ".wiff")[[1]][1]
}

# merge!

keys <- merge(sampkey, mousekey, by = "id")

# clean

toSave <- keys %>% 
  select(id, shortID, tissue.type, ion, treatmentID, Exercise, Weight)

write.table(keys, "allLipidKeys.txt", row.names = FALSE)

## Finish cleaning pos and neg data

# read without row1 because that is MetaboAnalyst "grouping" variable
pos <- read.csv("F:/Metabolomics/Rogers_TripleTOF/POS_MetaboAnalyst/data_normalized.csv")[-1,]
neg <- read.csv("F:/Metabolomics/Rogers_TripleTOF/NEG_MetaboAnalyst/data_normalized.csv")[-1,]

key <- read.table("allLipidKeys.txt", header = TRUE)

posv <- pos %>% 
  gather(key = longID, value = "area", 2:ncol(pos)) %>% 
  rename(lipid = X)
# add shortID column for grouping variable to merge
for(i in 1:nrow(posv)) {
  posv$shortID[i] <- strsplit(posv$longID[i], ".D")[[1]][1]
}


negv <- neg %>% 
  gather(key = longID, value = "area", 2:ncol(neg)) %>% 
  rename(lipid = X)
# add shortID column for grouping variable to merge
for(i in 1:nrow(negv)) {
  negv$shortID[i] <- strsplit(negv$longID[i], ".D")[[1]][1]
}

# merge all with the key
all <- rbind(negv, posv) %>% 
  merge(key, by = "shortID") %>% 
  # have lots of columns; don't need to keep al
  select(lipid, SampleName, tissue.type, ion, mouseID, treatmentID, Exercise, Weight, area)

# write this to file
write.table(all, "allLipidsCleaned.txt", row.names = FALSE)
