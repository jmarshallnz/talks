library(tidyverse)
library(islandR)

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

at_i <- attribution(ST ~ 1, st_i, data=all %>% filter(Source == "Human"))
at_d <- attribution(ST ~ 1, st_d, data=all %>% filter(Source == "Human"))

bind_rows(predict(at_i, FUN=identity) %>% mutate(Model = "Island"),
          predict(at_d, FUN=identity) %>% mutate(Model = "Dirichlet")) %>%
  group_by(Source, Model) %>% summarise(m = median(p), lci = quantile(p, 0.025),
                                        uci = quantile(p, 0.975), ucl = quantile(p, 0.75),
                                        lcl = quantile(p, 0.25)) %>% ungroup() -> foo

write_csv(foo, "data/overall_attribution.csv")
ggplot(foo, aes(x=Source, col=Model)) + geom_boxplot(aes(lower = lcl, upper = ucl,
                                              middle = m,
                                              ymin = lci, ymax = uci), stat="identity")

library(lubridate)
seasons <- tibble::tribble(~Month, ~Season,
                           3, "Autumn",
                           6, "Winter",
                           9, "Spring",
                           12, "Summer") %>%
  mutate(Season = fct_inorder(Season))

season_human <- all %>% filter(Source == "Human") %>% mutate(SampledDate = as.Date(SampledDate)) %>%
  mutate(Month = month(floor_date(SampledDate, "season"))) %>% left_join(seasons, by='Month') %>%
  filter(Year >= 2008)

at_is <- attribution(ST ~ Season*UR2006, st_i, data=season_human)
at_ds <- attribution(ST ~ Season*UR2006, st_d, data=season_human)

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
