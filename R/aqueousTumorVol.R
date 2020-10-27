
# ---- read data ----

require(tidyverse)

# define standard error function
se <- function(x) sqrt(var(x)/length(x))

dat <- read.table("https://raw.githubusercontent.com/EmilyB17/mice-metab/master/data/cleanedNamesAqueous.txt", sep = "\t", header = TRUE) %>% 
  # hyphenate instead of _ for treatment ID
  mutate(treatmentID = str_replace_all(treatmentID, "_", "-"))

# read tumor measurements
tum <- read.table("https://github.com/EmilyB17/mice-metab/raw/master/data/tumor-measurements.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# set ggplot theme 
theme_set(theme_bw())

# get colors
source("./R/RColorBrewer.R")

## wrangle to merge dat and tum
both <- dat %>% 
  mutate(newid = str_extract(id, "15-02-M(\\d){1,2}"),
         # add Time
         Time = case_when(
           str_detect(id, "_plasmaD7") ~ "D7",
           str_detect(id, "_plasmaD21") ~ "D21",
           str_detect(id, "_plasmaD35") ~ "D35",
           # for this, we need "tumor" to be "d35"
           str_detect(id, "_tumor") ~ "D35"
         )) %>% 
  left_join(tum, by = c("newid" = "Mouse", "Time"))

## ---- tumor regressions ----

# get tumor df
df <- both %>% filter(str_detect(id, "_tumor")) 

metab <- unique(both$Metabolite)
pvals <- data.frame()

for(i in 1:length(metab)) {
  
  mod <- summary(lm(area ~ cm3, data = filter(df, Metabolite == metab[i])))
  
  # get p values
  p <- data.frame(Metabolite = metab[i],
                  raw.p = round(mod$coefficients[2, 4], 3),
                  r2 = round(mod$r.squared, 3),
                  adjp = round(p.adjust(mod$coefficients[2, 4], method = "fdr"), 3),
                  slope = round(mod$coefficients[2,1], 3))
  
  # store p vals
  pvals <- rbind(pvals, p)
}

# get significant p values
sigs <- pvals %>% 
  filter(adjp < 0.05) %>% 
  # get column for positive and negative slopes
  mutate(sign = case_when(
    slope > 0 ~ "pos",
    slope < 0 ~ "neg"
  ))

# save table of results
#write.table(sigs, "./data/tumorvol-regression-tumor.txt", sep = "\t", row.names = FALSE)

# plot all
pdf <- df %>% 
  semi_join(sigs, by = "Metabolite")

ggplot(data = pdf, aes(x = area, y = cm3)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Metabolite, scales = "free") +
  labs(x = "Relative concentration", y = "Tumor volume")

# save 
#ggsave(filename = "./data/plots/tumorreg-tumor.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")

## ---- plasma regressions ----

# get plasma
df <- both %>% filter(str_detect(id, "_plasma")) %>% 
  # remove day7 since most volumes are 0
  filter(Time != "D7")

metab <- unique(df$Metabolite)
pvals <- data.frame()

for(i in 1:length(metab)) {
  
  # regression (this considers BOTH day 21 and day 35 in the same data)
  mod <- summary(lm(area ~ cm3, data = filter(df, Metabolite == metab[i])))
  
  # get p values
  p <- data.frame(Metabolite = metab[i],
                  raw.p = round(mod$coefficients[2, 4], 3),
                  r2 = round(mod$r.squared, 3),
                  adjp = round(p.adjust(mod$coefficients[2, 4], method = "fdr"), 3),
                  slope = round(mod$coefficients[2,1], 3))
  
  # store p vals
  pvals <- rbind(pvals, p)
}

# get significant p values
sigs <- pvals %>% 
  filter(adjp < 0.05) %>% 
  # get column for positive and negative slopes
  mutate(sign = case_when(
    slope > 0 ~ "pos",
    slope < 0 ~ "neg"
  ))


# save table of results
#write.table(sigs, "./data/tumorvol-regression-plasma-alltime.txt", sep = "\t", row.names = FALSE)

# plot all
pdf <- df %>% 
  semi_join(sigs, by = "Metabolite") %>% 
  mutate(time1 = case_when(
    Time %in% "D21" ~ "Day 21",
    Time %in% "D35" ~ "Day 35"
  ))

ggplot(data = pdf, aes(x = area, y = cm3)) +
  geom_point(aes(color = time1)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Metabolite, scales = "free") +
  labs(x = "Relative concentration", y = "Tumor volume") +
  scale_color_manual(values = timecols)

# save 
#ggsave(filename = "./data/plots/tumorreg-plasma.png", dpi = 600, plot = last_plot(), height = 4.35, width = 7.32, units = "in")


## ---- stepwise regression; tumor ----

require(MASS)

# get tumor df
df <- both %>% filter(str_detect(id, "_tumor")) %>% 
  pivot_wider(names_from = Metabolite, values_from = area) %>% 
  column_to_rownames(var = "newid") %>% 
  dplyr::select(-c(id, treatmentID, Exercise, Weight, Time))

dfmat <- as.matrix(df)
# use AIC as criterion
full.model <- lm(cm3 ~ ., data = df)
step.model <- stepAIC(full.model, direction = "both", trace = FALSE)
