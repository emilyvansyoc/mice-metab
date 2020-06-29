Metabolomics PermANOVA
================
Emily Bean
February 23, 2020

Overview
========

This script follows the same structure as `metabolmicsAnalysis.Rmd` but performs a permutational ANOVA on the entire metabolomics set. The data is scaled from 0 to 1, and then a Bray-Curtis distance matrix is created. A pairwise permutational ANOVA (PERMANOVA) on the data based on the grouping variable.

This analysis results in a pairwise comparison on the "structure" of the data as a "community".

Pairwise comparisons
====================

*All comparisons made for both tumor and plasma tissues*

1.  Plasma vs tumor
2.  Plasma D7 vs D21 vd D35
3.  4 treatment groups (2x2 factorial)
4.  Exercise vs sedentary
5.  Weight gain vs weight maintenance

Plasma vs Tumor
---------------

    ## $parent_call
    ## [1] "aqhdis ~ tissue.type , strata = Null"
    ## 
    ## $plasma_vs_tumor
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
    ## tissue.type  1    1.1981 1.19811  87.187 0.5123  0.001 ***
    ## Residuals   83    1.1406 0.01374         0.4877           
    ## Total       84    2.3387                 1.0000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

Plasma D7 vs D21 vs D35
-----------------------

    ## $parent_call
    ## [1] "plasmahdis ~ Time , strata = Null"
    ## 
    ## $`35_vs_7`
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
    ## Time       1   0.04900 0.048996  2.6517 0.05809  0.002 **
    ## Residuals 43   0.79453 0.018477         0.94191          
    ## Total     44   0.84352                  1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`35_vs_21`
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
    ## Time       1   0.03545 0.035448   2.128 0.04934  0.009 **
    ## Residuals 41   0.68298 0.016658         0.95066          
    ## Total     42   0.71843                  1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`7_vs_21`
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## Time       1   0.01921 0.019206   1.062 0.02719  0.334
    ## Residuals 38   0.68719 0.018084         0.97281       
    ## Total     39   0.70639                  1.00000       
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

Tumor
-----

### Four Treatment Groups (2x2 factorial)

    ## $parent_call
    ## [1] "tumorhdis ~ treatmentID , strata = Null"
    ## 
    ## $EX_AL_vs_EX_ER
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
    ## treatmentID  1   0.06600 0.065999  1.6734 0.19294  0.099 .
    ## Residuals    7   0.27608 0.039440         0.80706         
    ## Total        8   0.34208                  1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $EX_AL_vs_SED_ER
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
    ## treatmentID  1   0.05746 0.057455  1.6894 0.17436  0.074 .
    ## Residuals    8   0.27207 0.034009         0.82564         
    ## Total        9   0.32953                  1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $EX_AL_vs_SED_AL
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## treatmentID  1   0.01617 0.016171 0.25291 0.03064  0.859
    ## Residuals    8   0.51153 0.063941         0.96936       
    ## Total        9   0.52770                  1.00000       
    ## 
    ## $EX_ER_vs_SED_ER
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## treatmentID  1  0.011628 0.011628 0.72714 0.07475  0.647
    ## Residuals    9  0.143922 0.015991         0.92525       
    ## Total       10  0.155550                  1.00000       
    ## 
    ## $EX_ER_vs_SED_AL
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## treatmentID  1   0.03552 0.035518  0.8338 0.08479  0.622
    ## Residuals    9   0.38338 0.042597         0.91521       
    ## Total       10   0.41889                  1.00000       
    ## 
    ## $SED_ER_vs_SED_AL
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## treatmentID  1   0.03610 0.036102 0.95163 0.08689  0.504
    ## Residuals   10   0.37937 0.037937         0.91311       
    ## Total       11   0.41547                  1.00000       
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

### Exercise vs. Sedentary

    ## $parent_call
    ## [1] "tumorhdis ~ Exercise , strata = Null"
    ## 
    ## $EX_vs_SED
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
    ## Exercise   1   0.00657 0.006570 0.16477 0.0086  0.997
    ## Residuals 19   0.75755 0.039871         0.9914       
    ## Total     20   0.76412                  1.0000       
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

### Weight gain vs. weight maintenance

    ## $parent_call
    ## [1] "tumorhdis ~ Weight , strata = Null"
    ## 
    ## $AL_vs_ER
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
    ## Weight     1   0.08087 0.080872  2.2489 0.10584   0.04 *
    ## Residuals 19   0.68325 0.035960         0.89416         
    ## Total     20   0.76412                  1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

Plasma
------

### Exercise vs. Sedentary

    ## $parent_call
    ## [1] "plasmahdis ~ Exercise , strata = Null"
    ## 
    ## $EX_vs_SED
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
    ## Exercise   1    0.0190 0.019003  1.0393 0.01649  0.338
    ## Residuals 62    1.1336 0.018284         0.98351       
    ## Total     63    1.1526                  1.00000       
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"

### Weight gain vs. weight maintenance

    ## $parent_call
    ## [1] "plasmahdis ~ Weight , strata = Null"
    ## 
    ## $ER_vs_AL
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
    ## Weight     1   0.03281 0.032812  1.8167 0.02847  0.064 .
    ## Residuals 62   1.11980 0.018061         0.97153         
    ## Total     63   1.15262                  1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## attr(,"class")
    ## [1] "pwadstrata" "list"
