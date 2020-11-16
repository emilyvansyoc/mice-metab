


### Build errorplot with mean + standard error
p <- ggerrorplot(data = plotdat, x = "Time", y = "ab.change",
                 # mean and se
                 desc_stat = "mean_se", error.plot = "errorbar", width = 0.6,
                 # add dotplot and fill by Time
                 add = "dotplot", add.params = list(fill = "Time", size = 3),
                 # facet by Metabolite
                 facet.by = "Metabolite", scales = "free", ncol = 2,
                 # change axis titles
                 xlab = "Time", ylab = "\u0394 SED+AL") +
  # add mean as a line in the middle
  stat_summary(geom = "point", shape = 95, fun = "mean", col = "black", size = 10) +
  # color greyscale
  scale_fill_manual(values = timeGreys) +
  scale_y_continuous(expand = expansion(mult = 0, add = c(0, 0.6)))

ggpar(p, ggtheme = theme_pubr())+
  # make facet wrap title background white
  theme(strip.background = element_rect(
    fill="white", linetype=0
  ), 
  strip.text.x = element_text(size = 10, face = "bold"))


### REMOVE D21 and do T test between Day 7 and Day 35
datf <- abdat %>% filter(!Time == "Day 21") %>% filter(!treatmentID == "SED+AL")

# get unique metabolites and treatments
trt <- unique(datf$treatmentID)
mets <- unique(datf$Metabolite)

pvals <- data.frame()

for(j in 1:length(trt)) {
  
  for(i in 1:length(mets)) {
    
    # define model
    mod <- t.test(ab.change ~ Time, data = filter(datf, treatmentID == trt[j] & Metabolite == mets[i]))
    
    
    # get p vals
    ps <- data.frame(Metabolite = mets[i],
                     treatmentID = trt[j],
                     pval = round(mod$p.value, 3),
                     Tstat = round(mod$statistic, 3),
                     row.names = NULL)
    
    # concatenate
    pvals <- rbind(pvals, ps)
  }
}

# add FDR adjustment  
pvals$padj <- p.adjust(pvals$pval, method = "fdr")

# get sigs
sigs <- pvals %>% drop_na() %>% filter(padj < 0.05)
