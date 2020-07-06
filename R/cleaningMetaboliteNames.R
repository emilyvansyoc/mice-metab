
## Fix metabolite names based on KEGG classifications

require(tidyverse)
dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/allAqueousCleaned.txt", header = TRUE, stringsAsFactors = TRUE) %>% 
  mutate(tissue.type =
           case_when(Label %in% "tumor" ~ "tumor",
                     Label %in% "plasmaD7" ~ "plasma",
                     Label %in% "plasmaD21" ~ "plasma",
                     Label %in% "plasmaD35" ~ "plasma"))
source("./R/metabolomicsKEGG.Rmd")

# merge names - use "both" dataframe
df <- both %>% 
  mutate(Metabolite = case_when(
    is.na(Match) ~ metabolite,
    !is.na(Match) ~ Match
  )) %>% 
  # remove extra columns
  select(-c(metabolite, Match, HMDB, PubChem, KEGG, SMILES))

# write to file
write.table(df, file = "./data/cleanedNamesAqueous.txt", sep = "\t", row.names = FALSE)
