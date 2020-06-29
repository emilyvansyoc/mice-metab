

## script to load or install the normally required packages

## ---- CRAN packages ----

cran.packages <- c("ggplot2", 
                   "vegan", 
                   "devtools", 
                   "tidyr", 
                   "dplyr",
                   "emmeans")
not_installed <- cran.packages[!(cran.packages %in% installed.packages()[ , "Package"])] 

if(length(not_installed)) {
  
  # print warning message that dependencies install
  cat("\n installing packages:", not_installed, "\n")
  # use error handling to install CRAN packages
  try(install.packages(not_installed), silent = TRUE)
  
 
  
}

reqs <- paste0("require(", cran.packages, ")")

tryCatch(c(require(ggplot2),
           require(vegan),
           require(devtools),
           require(tidyr),
           require(dplyr),
           require(emmeans)),
         
         error = function(c) "Packages are not installed or loaded correctly")

### ---- Bioconductor ----

bioc.packages <- c("ggordplots",
                   "pairwiseAdonis",
                   "phyloseq",
                   "dada2")

not_installed <- bioc.packages[!(bioc.packages %in% installed.packages()[ , "Package"])] 


if(length(not_installed)) {
  
  # print warning message that dependencies install
  cat("\n installing packages:", not_installed, "\n")
 
  # use error handling to install github/bioconductor packages
  require("devtools")
  try(install_github("jfq3/ggordiplots"), silent = TRUE)
  try(install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"), silent = TRUE)
  
  
}


tryCatch(c(require(phyloseq),
           require(dada2),
           require(ggordiplots),
           require(pairwiseAdonis)),
           
         
         error = function(c) "Packages are not installed or loaded correctly")
