library(tidyverse)
library(islandR)
library(lubridate)

# Fit the sequence type distribution using the island model
num_samples <- 100

set.seed(2)

all = read.csv("data/attribution_data.csv")

st_i = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = all,
              method = "island",
              sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
              samples=num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = all,
              method = "dirichlet",
              samples = num_samples)

seasons <- tibble::tribble(~Month, ~Season,
                           3, "Autumn",
                           6, "Winter",
                           9, "Spring",
                           12, "Summer") %>%
  mutate(Season = fct_inorder(Season))

humans <- all %>% filter(Source == "Human") %>% mutate(SampledDate = as.Date(SampledDate)) %>%
  mutate(Month = month(floor_date(SampledDate, "season"))) %>% left_join(seasons, by='Month') %>%
  filter(Year >= 2008) %>% mutate(UR2006_fact = factor(UR2006_num))

at_i <- attribution(ST ~ 1, st_i, data=humans)
at_d <- attribution(ST ~ 1, st_d, data=humans)

bind_rows(predict(at_i, FUN=identity) %>% mutate(Model = "Island"),
          predict(at_d, FUN=identity) %>% mutate(Model = "Dirichlet")) %>%
  group_by(Source, Model) %>% summarise(m = median(p), lci = quantile(p, 0.025),
                                        uci = quantile(p, 0.975), ucl = quantile(p, 0.75),
                                        lcl = quantile(p, 0.25)) %>% ungroup() -> foo

write_csv(foo, "data/overall_attribution.csv")
ggplot(foo, aes(x=Source, col=Model)) + geom_boxplot(aes(lower = lcl, upper = ucl,
                                              middle = m,
                                              ymin = lci, ymax = uci), stat="identity")

at_is <- attribution(ST ~ Season*UR2006, st_i, data=humans)
at_ds <- attribution(ST ~ Season*UR2006, st_d, data=humans)

bind_rows(predict(at_is, FUN=identity) %>% mutate(Model = "Island"),
       predict(at_ds, FUN=identity) %>% mutate(Model = "Dirichlet")) %>%
  extract(X, into="Season", regex="Season([^:]*)=1", remove=FALSE) %>%
  extract(X, into="UR2006", regex="UR2006([^:]*)=1") %>%
  replace_na(list(Season = "Autumn", UR2006 = "Rural")) %>%
  group_by(Source, UR2006, Season, Model) %>% summarise(m = median(p), lci = quantile(p, 0.025),
                                        uci = quantile(p, 0.975), ucl = quantile(p, 0.75),
                                        lcl = quantile(p, 0.25)) %>% ungroup() -> foo

write_csv(foo, "data/season_attribution.csv")

ggplot(foo, aes(x=Source, col=UR2006)) + geom_boxplot(aes(lower = lcl, upper = ucl,
                                                         middle = m,
                                                         ymin = lci, ymax = uci), stat="identity") +
  facet_grid(Model~Season)

at_if <- attribution(ST ~ UR2006_fact, st_i, data=humans)
at_df <- attribution(ST ~ UR2006_fact, st_d, data=humans)

at_il <- attribution(ST ~ UR2006_num, st_i, data=humans)
at_dl <- attribution(ST ~ UR2006_num, st_d, data=humans)

df_if <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR2006_num=-3)) %>% mutate(GenotypeModel = 'Island', AttributionModel='Categorical')
df_il <- predict(at_il, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_df <- predict(at_df, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR2006_num=-3)) %>% mutate(GenotypeModel = 'Dirichlet', AttributionModel='Categorical')
df_dl <- predict(at_dl, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

bind_rows(df_if, df_il, df_df, df_dl) %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num) %>%
  summarise(m = median(p), lci = quantile(p, 0.025),
                                                        uci = quantile(p, 0.975), ucl = quantile(p, 0.75),
                                                        lcl = quantile(p, 0.25)) %>% ungroup() -> foo

write_csv(foo, "data/cat_lin_compare.csv")
ggplot(foo) +
  geom_ribbon(aes(x=UR2006_num, ymin=lc, ymax=uc, fill=Source), alpha=0.3) +
  geom_ribbon(aes(x=UR2006_num, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR2006_num, y=m, col=Source), lwd=1) +
  facet_grid(GenotypeModel~AttributionModel) +
  scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
  scale_y_continuous(name="Percentage of cases", labels=scales::percent_format(), limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
  theme(axis.text.x = element_text(hjust=c(-0.1,1.1)),
        axis.text.y = element_text(vjust=c(-0.1,rep(0.5, 3), 1.1)))

