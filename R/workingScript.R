
# separating positive and negative ions

pos <- lip %>% filter(ion == "POS") %>% 
  mutate(lipname = sapply(str_split(lipid, ";"), `[`, 1))

neg <- lip %>% filter(ion == "NEG") %>% 
  mutate(lipname = sapply(str_split(lipid, ";"), `[`, 1))

unique(pos$lipid)
unique(neg$lipid)

posn <- paste(unique(pos$lipname), sep = "|")
negn <- unique(neg$lipname)

# there are 205 that are in both positive and negative runs
match <- posn[posn %in% negn]


### interaction plots

# get metabolites with significant interaction term
sigs <- data.frame(plasmamod[1]) %>% 
  filter(pvalEx_Wt < 0.05)

sigdf <- plasma %>% 
  semi_join(sigs, by = "metabolite") %>% 
  group_by(metabolite, Exercise, Weight) %>% 
  summarize(mean = mean(area),
            sd = sd(area),
            cilow = mean - 1.96*sd,
            cihi = mean + 1.96*sd) %>% 
  ungroup()

# make a plot
ggplot(data = sigdf, aes(x = Exercise, y = mean, group = Weight, color = Weight)) +
  geom_point() +
  geom_line() +
  geom_errorbar(data = sigdf, aes(ymin = mean + sd, ymax = mean - sd), width = .1) +
  facet_wrap(~metabolite) +
  ggtitle("Significant Interactions w/o Time") +
  labs(x = "Exercise", y = "Mean area")

## interactions with Time
sigs <- data.frame(plasmamod[1]) %>% 
  filter(pvalEx_Wt_Time < 0.05)

sigdf <- plasma %>% 
  semi_join(sigs, by = "metabolite") %>% 
  group_by(metabolite, Exercise, Weight, Time) %>% 
  summarize(mean = mean(area),
            sd = sd(area),
            cilow = mean - 1.96*sd,
            cihi = mean + 1.96*sd) %>% 
  ungroup()

# make plots
ggplot(data = sigdf) +
  geom_point(aes(x = Time, y = mean, group = Exercise, color = Exercise)) +
  geom_line(aes(x = Time, y = mean, group = Exercise, color = Exercise)) +
  geom_point(aes(x = Time, y = mean, group = Weight, color = Weight)) +
  geom_line(aes(x = Time, y = mean, group = Weight, color = Weight))
  #geom_line() +
  #geom_errorbar(data = sigdf, aes(ymin = mean + sd, ymax = mean - sd), width = .1) +
  facet_wrap(~metabolite) +
  ggtitle("Significant Interactions") +
  labs(x = "Exercise", y = "Mean area")


ggplot(data = sigdf, aes(x = Time, y = mean, group = treatmentID, color = treatmentID)) +
  geom_point(aes(shape = Exercise)) +
  geom_line(aes(linetype = Weight)) +
  facet_wrap(~metabolite) +
  ggtitle("Significant Interactions")+
  labs(x = "Exercise", y = "Mean area")


sigdf <- plasma %>% 
  semi_join(sigs, by = "metabolite") %>% 
  group_by(metabolite, treatmentID, Time, Exercise, Weight) %>% 
  summarize(mean = mean(area),
            sd = sd(area),
            cilow = mean - 1.96*sd,
            cihi = mean + 1.96*sd) %>% 
  ungroup()


## get p values from tukey instead

sig <- data.frame(tumormod[2])
