# define a function to perform two-way ANOVA
myTwoWayAnova <- function(df, metabList) {
  
  
  # make output dataframes
  allpvals <- data.frame()
  posthoc_sigs <- data.frame()
  
  # loop through each metabolite
  for(i in 1:length(metabList)) {
    
    # perform ANOVA on the iterated metabolite
    mod <- aov(area ~ Exercise*Weight, data = filter(df, metabolite == metabList[i]))
    
    # collect p values and F statistics
    pvals <- data.frame(metabolite = metabList[i],
                        fStatEx =  round(summary(mod)[[1]][["F value"]][1], 3),
                        pvalEx = round(summary(mod)[[1]][["Pr(>F)"]][1], 3),
                        fStatWt = round(summary(mod)[[1]][["F value"]][2], 3),
                        pvalWt = round(summary(mod)[[1]][["Pr(>F)"]][2], 3),
                        fStatInt = round(summary(mod)[[1]][["F value"]][3], 3),
                        pvalInt = round(summary(mod)[[1]][["Pr(>F)"]][3], 3))
    
    # save to an outDF
    allpvals <- rbind(allpvals, pvals)
    
    # perform Tukey post-hoc
    tukey <- TukeyHSD(mod)
    
    # collect p vals
    posthoc <- as.data.frame(rbind(as.data.frame(tukey[1]$Exercise) %>% 
                                     mutate(comparison = "Exercise",
                                            contrast = rownames(tukey[1]$Exercise)),
                                   as.data.frame(tukey[2]$Weight) %>% 
                                     mutate(comparison = "Weight",
                                            contrast = rownames(tukey[2]$Weight)),
                                   as.data.frame(tukey[3]$`Exercise:Weight`) %>% 
                                     mutate(comparison = "Interaction",
                                            contrast = rownames(tukey[3]$`Exercise:Weight`))))
    
    
    # if there are any significant p vals, add to outDF
    if(any(posthoc$`p adj` < 0.05)) {
      
      out <- posthoc %>% filter(`p adj` < 0.05)
      out$metabolite <- metabList[i]
      
      posthoc_sigs <- rbind(posthoc_sigs, out)
      
    }
    
    
  }
  
  # return the output dataframes
  return(list(allpvals, posthoc_sigs))
}
