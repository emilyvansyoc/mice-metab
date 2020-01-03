## LIPIDS: Data Wrangling
require(dplyr)
require(tidyr)
require(stringr)

sampkey <- read.csv("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/lipidSampleKeys.csv", stringsAsFactors = FALSE)[-1,] # remove Blank row

mousekey <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/mouseKeys.txt", header = TRUE, stringsAsFactors = FALSE)


## wrangle to add the proper identifying variables

# add mouse ID, exercise and weight IDs
for(i in 1:nrow(mousekey)) {
  mousekey[i, "id"] <- strsplit(mousekey$mouseID[i], "-")[[1]][3]
  mousekey[i, "Exercise"] <- strsplit(mousekey$treatmentID[i], "_")[[1]][1]
  mousekey[i, "Weight"] <- strsplit(mousekey$treatmentID[i], "_")[[1]][2]
}

for(i in 1:nrow(sampkey)) {
  
  # if the SampleName has "_" in it, it is a plasma sample (ex D35_M1)
  if(grepl("_", sampkey$SampleName[i]) == TRUE) {
    
    # separate out the mouse ID
    sampkey[i, "id"] <- strsplit(sampkey$SampleName[i], "_")[[1]][2]
    sampkey[i, "tissue.type"] <- "plasma"
   
    # if the SampleName is only one variable, it is a tumor sample (ex M7)
  } else {
    
    sampkey[i, "id"] <- sampkey$SampleName[i]
    sampkey[i, "tissue.type"] <- "tumor"
    
  }
  
  # key out positive and negative
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

write.table(keys, "./data/allLipidKeys.txt", row.names = FALSE)

## Finish cleaning pos and neg data

# read without row1 because that is MetaboAnalyst "grouping" variable
pos <- read.csv("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/data_normalizedPOS.csv", stringsAsFactors = FALSE)[-1,]
neg <- read.csv("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/data_normalizedNEG.csv", stringsAsFactors = FALSE)[-1,]

##key <- read.table("allLipidKeys.txt", header = TRUE)

# make the lipid spreadsheets vertical for merging

# collect the column names minus the first one ("X", lipid identifier)
ids <- colnames(pos)[-1]

posv <- pos %>% 
  gather(key = longID, value = "log.area", ids) %>% 
  rename(lipid = X)
# add shortID column for grouping variable to merge
for(i in 1:nrow(posv)) {
  posv$shortID[i] <- strsplit(posv$longID[i], "\\.")[[1]][1]
  posv$SampleName[i] <- strsplit(posv$longID[i], "\\.")[[1]][2]
}

negids <- colnames(neg)[-1]
negv <- neg %>% 
  gather(key = longID, value = "log.area", negids) %>% 
  rename(lipid = X)
# add shortID column for grouping variable to merge
for(i in 1:nrow(negv)) {
  negv$shortID[i] <- strsplit(negv$longID[i], "\\.")[[1]][1]
  negv$SampleName[i] <- strsplit(negv$longID[i], "\\.")[[1]][2]
}

# merge all with the key
all <- rbind(negv, posv) %>% 
  merge(keys, by = c("SampleName", "shortID"))  %>% 
  # have lots of columns; don't need to keep all
  select(lipid, SampleName, tissue.type, ion, mouseID, treatmentID, Exercise, Weight, log.area)

# write this to file
write.table(all, "./data/allLipidsCleaned.txt", row.names = FALSE)
