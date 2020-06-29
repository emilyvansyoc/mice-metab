


### CUSTOM FUNCTION TO RUN COMMUNITY ANALYSIS

# input: Two dataframes, one with abundance data and one with sample variables -- these 
# MUST BE ORDERED THE SAME WAY

# options: select the variable to test with ADONIS and the number of NMDS dimensions

# the function performs normalization to relative abundance, an adonis statistical test,
# and a pairwise adonis test (if the first adonis test is significant) with bonferroni correction. 
# Then, an ordination is created and a stressplot is printed. A basic NMDS is created
# using two packages: ggplot2 and ggordiplots

myNMDS <- function(abundDF, # dataframe with only abundance data
                   sampDF, # dataframe with sample information and test variables
                   indVar, # variable to test
                   ordK = 2 # number of dimensions for ordination
) {
  
  ## ---- error handling for requiring packages ----
  
  need.packages <- c("pairwiseAdonis", "vegan", "ggplot2", "ggordiplots", "devtools")
  not_installed <- need.packages[!(need.packages %in% installed.packages()[ , "Package"])] 
  if(length(not_installed)) {
    
    # print warning message that dependencies install
    cat("\n installing packages:", not_installed, "\n")
    # use error handling to install CRAN packages
    try(install.packages(not_installed), silent = TRUE)
    
    # use error handling to install github packages
    require("devtools")
    try(install_github("jfq3/ggordiplots"), silent = TRUE)
    try(install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"), silent = TRUE)
    
  }
  
  tryCatch(c(require(pairwiseAdonis),
             require(vegan),
             require(ggplot2),
             require(ggordiplots)),
           
           error = function(c) "Packages are not installed or loaded correctly")
  
  ## ---- stat testing ----
  
  ## make relative abundance
  # use decostand from vegan
  relabun <- decostand(abundDF, method = "range", MARGIN = 2)
  
  ## make distance matrix
  # use vegdist from vegan
  dis <- vegdist(relabun, method = "euclidean")
  
  ## adonis test
  # create formula
  f <- as.formula(paste("dis ~", indVar))
  print(adonis(f, data = sampDF)) # this will print to the console
  
  # if significant: run pairwise adonis
  a <- adonis(f, data = sampDF) 
  if(a$aov.tab$`Pr(>F)`[1] < 0.05) {
    
    cat("\n printing pairwise adonis test with bonferroni correction \n")
    print(pairwise.adonis2(f, data = sampDF, p.adjust.m = "bon", perm = 1000))
    
  } else {
    
    warning("\n adonis test: no variables were significant (p < 0.05) \n")
  }
  
  ## make ordination
  ord <- metaMDS(relabun, k = ordK)
  # print stress plot
  stressplot(ord, dis, main = paste0("Stressplot of ordination fit, dimensions = ", ordK))
  
  ### ---- create plot ----
  
  # get points from ordination
  sampDF$NMDS1 <- ord$points[,1]
  sampDF$NMDS2 <- ord$points[,2]
  
  # save gg_ordiplot object to get ellipse values
  plot <- gg_ordiplot(ord, groups = sampDF[, indVar], label = FALSE, plot = FALSE)
  # get ellipse coordinates
  df_ell <- plot$df_ellipse
  # get label coordinates for ellipse centers
  NMDS.mean <- plot$df_mean.ord
  # pull NMDS coordinates
  ord.data <- plot$df_ord
  # create in ggplot2
  gp <- ggplot(data = sampDF, aes(x = NMDS1, y = NMDS2)) +
    geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
    geom_point(aes(x = NMDS1, y = NMDS2, color = sampDF[, indVar]), size = 1) +
    annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group) +
    ggtitle("NMDS Ordination") +
    theme_bw()
  
  print(gp)
  
  
}



