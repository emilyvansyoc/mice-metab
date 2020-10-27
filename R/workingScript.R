# get list of metabolites in table format for Supp Materials

mets <- unique(dat$Metabolite)
write.table(mets, file = "./data/aqueousMetaboliteList.txt", sep = "\t", row.names = FALSE)
