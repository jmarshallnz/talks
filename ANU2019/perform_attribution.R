library(islandR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
# filter out those where we don't have a location
attr <- read.csv("data/extract_attribution.csv")

# NOTE: Bug in fct_collapse in forcats vs 0.4.0. Fixed with the Other line
dataset = attr %>% mutate(Source = fct_collapse(Source, Human = c("Human"), Ruminants = c("Cattle", "Sheep"),
                                      Poultry = c("Supplier A", "Supplier B", "Supplier Other"),
                                      Water = "Environmental water",
                                      Other = c("Cat_dog_pet", "Duck_poultry",
                                                "Pig", "Spent_hen", "Turkey",
                                                "Water_bird_wild", "Wild_bird_other"))) %>%
  filter(Source != "Human" | !is.na(UR2006_num)) %>%
  mutate(Year = year(SampledDate)) %>%
  mutate(Intervention = ifelse(Year < 2008, "Before", "After")) %>%
  mutate(UR2006_fact = factor(UR2006_num)) %>%
  mutate(Age = cut(Age, breaks = c(-Inf, 5, Inf), right = FALSE))

dataset %>% count(Source)

write_csv(dataset, "data/attribution_data.csv")

num_samples = 100

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "island",
              sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
              samples = num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "dirichlet",
              samples = num_samples)

# do the attribution
humans <- dataset %>% filter(Source == "Human")

# Dirichlet, island numeric vs factorial
at1 <- attribution(ST ~ UR2006_num, st_i, humans)
at2 <- attribution(ST ~ UR2006_num, st_d, humans)
at3 <- attribution(ST ~ UR2006_fact, st_i, humans)
at4 <- attribution(ST ~ UR2006_fact, st_d, humans)

# Dirichlet, island, numeric by intervention
at5 <- attribution(ST ~ UR2006_num*Intervention, st_i, humans)
at6 <- attribution(ST ~ UR2006_num*Intervention, st_d, humans)

# Dirichlet, island, numeric by age post 2008 only
post2008 <- humans %>% filter(Intervention == "After") %>%
  mutate(Age = cut(Age, breaks = c(-Inf, 5, Inf), right = FALSE))
at7 <- attribution(ST ~ UR2006_num*Age, st_i, post2008)
at8 <- attribution(ST ~ UR2006_num*Age, st_d, post2008)

save(st_i, st_d, at1, at2, at3, at4, at5, at6, at7, at8, file="attribution.Rdata")


df_il <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_dl <- predict(at_id, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_il, df_dl)

save(list=c('df', 'at_if', 'at_id'), file='attribution_fits_intervention_prec1.Rdata')

at_if <- attribution(ST ~ UR2006_num*Intervention, st_i, humans, priors=list(theta_mean=0, theta_prec=10))
at_id <- attribution(ST ~ UR2006_num*Intervention, st_d, humans, priors=list(theta_mean=0, theta_prec=10))

df_il <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_dl <- predict(at_id, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_il, df_dl)

save(list=c('df', 'at_if', 'at_id'), file='attribution_fits_intervention_prec10.Rdata')

load('attribution_fits_intervention.Rdata')

plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num, Intervention) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    lc = quantile(p, 0.25),
    uc = quantile(p, 0.75),
    ui = quantile(p, 0.9)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))),
         Intervention = factor(Intervention, levels=c("Before", "After"), labels=c("2005-2007", "2008-2014")))

set.seed(3)
iters <- sample(unique(df_island$identity), 20)
df_island <- df %>% filter(GenotypeModel == "Island", identity %in% iters)
write.csv(df_island, "iters_intervention.csv", row.names=FALSE)

# plot
pdf("fig_intervention.pdf", width=7, height=5)
ggplot(plot_df) +
  geom_ribbon(aes(x=UR2006_num, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR2006_num, y=m, col=Source), lwd=1) +
  facet_grid(GenotypeModel~Intervention) +
  scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
  scale_y_continuous(name="Percentage of cases", labels=scales::percent_format(), limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  theme_bw() +
  theme(text = element_text(family="Times"),
        legend.position = c(0.99,0.89),
        legend.justification = "right",
        legend.margin=margin(0,0,0,0),
        legend.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(hjust=c(-0.1,1.1)),
        axis.text.y = element_text(vjust=c(-0.1,rep(0.5, 3), 1.1)))
dev.off()
