
# WORKING: PERMANOVA

require(vegan)
require(pairwiseAdonis)

# standardize data

test <- decostand(aqh[, 8:ncol(aqh)], method = "range", MARGIN = 2)

# make distance matrix
dis <- vegdist(test, method = "bray")

# perform adonis
mod <- adonis(dis ~ aqh$tissue.type)

# perform pairwsise adonis
mod <- pairwise.adonis2(dis ~ treatmentID, data = aqh, p.adjust.m = "bon", perm = 1000)


# read all aqueous data
aq <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/allAqueousCleaned.txt", header = TRUE, stringsAsFactors = TRUE) %>% 
  mutate(tissue.type =
           case_when(Label %in% "tumor" ~ "tumor",
                     Label %in% "plasmaD7" ~ "plasma",
                     Label %in% "plasmaD21" ~ "plasma",
                     Label %in% "plasmaD35" ~ "plasma"))

## Make all data horizontal for PCAs
aqh <- aq %>% 
  spread(key = metabolite, value = area)

# standardize the data
aqhstand <- decostand(aqh[ , 8:ncol(aqh)], method = "range", MARGIN = 2)

# make a distance matrix
aqhdis <- vegdist(aqhstand, method = "bray")

## plasma only
plasmah <- aqh %>% 
  filter(tissue.type == "plasma") %>% 
  mutate(Time = factor(case_when(
    Label %in% "plasmaD7" ~ "7",
    Label %in% "plasmaD21" ~ "21",
    Label %in% "plasmaD35" ~ "35"
  ), ordered = TRUE, levels = c("7", "21", "35")))

# this data has an extra column; remove for standardization
plasmst <- plasmah %>% select(-c(id, Label, mouseID, treatmentID, Exercise, Weight, tissue.type, Time))

# standardize the data
plasmahstand <- decostand(plasmst, method = "range", MARGIN = 2)

# make a distance matrix
plasmahdis <- vegdist(plasmahstand, method = "bray")

## tumor only
tumorh <- aqh %>% 
  filter(tissue.type == "tumor")

# standardize the data
tumorhstand <- decostand(tumorh[ , 8:ncol(tumorh)], method = "range", MARGIN = 2)

# make a distance matrix
tumorhdis <- vegdist(tumorhstand, method = "bray")


## WORKING

pairwise.adonis2(plasmahdis ~ Time, data = plasmah)
