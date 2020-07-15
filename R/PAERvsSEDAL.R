
## Treatment versus Control analysis

# because PA+ER was found to have protective effects, what is changing in comparison to SED+AL?

# get data
require(tidyverse)
source("./R/RColorBrewer.R")
# treatments
trts <- c("EX_ER", "SED_AL")
dat <-  read.table("./data/aqueousCleanedNotNormalized.txt", sep = "", header = TRUE) %>% 
  filter(treatmentID %in% trts)

tumor <- dat %>% filter(Label %in% "tumor") %>% 
  select(-c(Exercise, Weight, id, Label))

plasma <- dat %>% filter(Label != "tumor") %>% 
  select(-c(Exercise, Weight, id))

## ---- T tests ----

# first, perform iterative t tests with adjusted p value to get any metabolites w/sig diffs
mets <- unique(tumor$metabolite)

sigs <- list()
nosig <- list()
for(i in 1:length(mets)) {
  
  # t test
  mod <- t.test(area ~ treatmentID, data = filter(tumor, metabolite == mets[i]))
  # adjust p value
  pad <- p.adjust(mod$p.value, method = "fdr", n = nrow(filter(tumor, metabolite == mets[i])))
  
  # get list of mets with significant p values
  if(pad < 0.05) {
    sigs <- c(sigs, mets[i])
  } else{
    nosig <- c(nosig, mets[i])
  }

}

# Interestingly, there are no significant metabolites 