

# get only plasma and calculate Time column
plasma <- dat %>% 
  filter(str_detect(id, "_plasma")) %>% 
  # make Time column
  mutate(Time = case_when(
    str_detect(id, "_plasmaD7") ~ "Day 7",
    str_detect(id, "_plasmaD21") ~ "Day 21",
    str_detect(id, "_plasmaD35") ~ "Day 35"
  ))

# remove Day 21
tdat <- plasma %>% filter(!Time == "Day 21")

# get unique metabolites and treatments
trt <- unique(tdat$treatmentID)
mets <- unique(tdat$Metabolite)
# L-Lactic acid is not found in the plasma; remove 
mets <- mets[-3]

pvals <- data.frame()

for(j in 1:length(trt)) {
  
  for(i in 1:length(mets)) {
    
    # define model
    mod <- t.test(area ~ Time, data = filter(tdat, treatmentID == trt[j] & Metabolite == mets[i]))
    
    # get p vals
    ps <- data.frame(Metabolite = mets[i],
                     treatmentID = trt[j],
                     pval = round(mod$p.value, 3),
                     tstat = round(mod$statistic, 3),
                     row.names = NULL)
    
    # concatenate
    pvals <- rbind(pvals, ps)
  }
}

# get significant 
sigs <- pvals %>% 
  drop_na() %>% 
  filter(pval < 0.05)

# get only PA+ER
paer <- sigs %>% filter(treatmentID == "PA+ER")

## PLOTS

# filter data
plotdat <- tdat %>% 
  semi_join(paer, by = "Metabolite") %>% 
  filter(treatmentID == "PA+ER") %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("Day 7", "Day 35")))

fig3 <- function(MetaboliteName) {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "Time", y = "ab.change", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "Time", width = 0.4, binwidth = 0.2),
                 # change axis titles
                 xlab = "Time", ylab = "\u0394 SED+AL",
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#525252")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              ggtheme = theme_pubr()) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))+
    rremove("x.text") +
    rremove("xlab") 
  
  return(p1)
  
}

# define a similar function to include labels for the bottom plots
fig3Labs <- function(MetaboliteName) {
  
  # define dataframe
  dat <- plotdat %>% filter(Metabolite == as.character(MetaboliteName))
  
  # define plot
  p <- ggbarplot(data = dat, x = "Time", y = "area", fill = "#D9D9D9",
                 #facet.by = "Metabolite", scales = "free",
                 add = c("mean_se", "dotplot"), add.params = list(fill = "Time", width = 0.4, binwidth = 0.2),
                 # change axis titles
                 xlab = "Time", ylab = "Relative Concentration",
                 # add title
                 title = as.character(MetaboliteName)) +
    scale_fill_manual(values = c("#FFFFFF", "#525252")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 0.5))) 
  
  p1 <- ggpar(p, legend = "none",
              ggtheme = theme_pubr()) +
    # center the plot title
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
  
}

met <- unique(plotdat$Metabolite)

p1 <- fig3Labs(met[1])
p2 <- fig3Labs(met[2])
p3 <- fig3Labs(met[3])
p4 <- fig3Labs(met[4])

g1 <- ggarrange(p1, 
                p2, 
                p3,
                p4,
                ncol = 2, nrow = 2,
                labels = c("A", "B", "C", "D"))
g1
