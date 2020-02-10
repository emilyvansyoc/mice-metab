## lipid data; aggregate by groups for pairwise comparisons
# need non-normalized data, then aggregate, then upload to MetaboAnalyst to transform & normalize

## ---- read and clean ----

# read data without first row (double ID row... pointless)
key <- read.table("./data/allLipidKeys.txt", header = TRUE, sep = "\t")
neg <- read.csv("./data/data_originalNEG.csv",
                  stringsAsFactors = FALSE, header = TRUE)[-1,] 
pos <- read.csv("./data/data_originalPOS.csv",
                stringsAsFactors = FALSE, header = TRUE)[-1,]

## WRANGLE
# collect the column names minus the first one ("X", lipid identifier)
ids <- colnames(pos)[-1]

posv <- pos %>% 
  gather(key = longID, value = "area", ids) %>% 
  rename(lipid = X) 


# same procedure for negative ions
negids <- colnames(neg)[-1]
negv <- neg %>% 
  gather(key = longID, value = "area", negids) %>% 
  rename(lipid = X)

### ---- aggregating by group ----

# there are 1329 positive lipids
length(unique(posv$lipid))
# and 1105 negative lipids 
length(unique(negv$lipid))

# get unique names
posv$group <- sapply(strsplit(posv$lipid, " "), `[`, 1)
posWo <- posv %>% filter(group == "w/o")
posNoWo <- posv %>% filter(!group == "w/o")

# get sum of each
posSum <- posNoWo %>% 
  mutate(area = as.numeric(area)) %>% 
  group_by(longID, group) %>% 
  summarize(areasum = sum(area),
            rangeLow = range(area)[1],
            rangeHigh = range(area)[2],
            count = length(unique(lipid)),
            type = "known")

# get unique w/o names (extra space in strsplit) 
posWo$noWo <- sapply(str_remove(posWo$lipid, "w/o "), `[`, 1)
posWo$group <- sapply(strsplit(posWo$noWo, " "), `[`, 1)

posSumWo <- posWo %>% 
  mutate(area = as.numeric(area)) %>% 
  group_by(longID, group) %>% 
  summarize(areasum = sum(area),
            rangeLow = range(area)[1],
            rangeHigh = range(area)[2],
            count = length(unique(lipid)),
            type = "w/o")
posSumWo$newGroup <- paste("w/o ", posSumWo$group)
posSumMerge <- posSumWo %>% 
  mutate(group = newGroup,
         newGroup = NULL)

# merge and clean for MetaboAnalyst
allpos <- rbind(posSum, posSumMerge) %>% 
  select(longID, group, areasum) %>% 
  spread(key = longID, value = areasum)

# write to file
write.table(allpos, "./data/posOriginalAggregated.csv", sep = ",", row.names = FALSE)

            