# anovas and kegg visualizations

## replace _ with - for plots
assignv <- assignv %>% 
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

theme_set(theme_minimal())

## tumor - barplot of metabolic classes
ggplot(data = filter(assignv, tissue.type == "tumor"), 
       aes(x = Class, y = area, fill = treatmentID)) +
  geom_col() +
  scale_fill_manual(values = treatmentIntcols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Metabolic Class", y = "Area Under the Curve", fill = "Treatment")

#ggsave(filename = "./data/plots/barplot-tumor-treatment-class.tiff", plot = last_plot(), dpi = "print")

## plasma - barplot of metabolic classes
ggplot(data = filter(assignv, tissue.type == "plasma"), 
       aes(x = Class, y = area, fill = treatmentID)) +
  geom_col() +
  scale_fill_manual(values = treatmentIntcols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Metabolic Class", y = "Area Under the Curve", fill = "Treatment")

#ggsave(filename = "./data/plots/barplot-plasma-treatment-class.tiff", plot = last_plot(), dpi = "print")

## test plot
test <- assignv %>% 
  # make a smaller dataset
  slice_sample(n = 10)

ggplot(data = test, aes(x = Exercise, y = area, fill = Metabolite)) +
  geom_col()

# get significant data
se <- function(x) sqrt(var(x)/length(x))

df <- as.data.frame(tumormod[2])

test <- assignv  %>% 
  filter(tissue.type == "tumor") %>% 
  semi_join(df, by = c("Metabolite" = "metabolite")) %>% 
  group_by(Metabolite, treatmentID) %>% 
  summarize(meanarea = mean(area),
            searea = se(area))
  
ggplot(data = test, aes(x = treatmentID, y = meanarea, fill = Metabolite)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = meanarea - searea, ymax = meanarea + searea), 
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_manual(values = accentcols)
