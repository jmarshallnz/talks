library(islandR)
library(dplyr)
library(lubridate)
library(tidyr)

#db = extract_for_attribution()
db = read.csv("data/sa_report/extract_attribution.csv", stringsAsFactors = FALSE)

# use a source map to map sources to final sources for attribution
source_map = read.csv("data/sa_report/5_source.csv")

db = db %>% left_join(source_map, by="Source") %>%
            filter(!is.na(Label))

# create a year column, and filter out years
db = db %>% mutate(Year = year(SampledDate)) %>%
  filter(Year >= 2005, Year <= 2016) %>%
  mutate(Year = ifelse(Year < 2008, 'Before', 'After'))

# fit the model
st = st_fit(Label ~ ST,
            sequences= ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
            non_primary="Human",
            method="island",
            data=db,
            iters=20000)

db$Year = as.factor(db$Year)
mod = attribution(ST ~ Year, sampling_dist = st, data=db %>% filter(Source == "Human"))

counts <- db %>% filter(Source=='Human') %>% group_by(Year) %>% summarize(Count=n()) %>%
  mutate(Years = ifelse(Year == "Before", 2.75, 9)) %>%
  mutate(Count = Count / Years) %>% select(Year, Count)

pred = predict(mod, FUN=identity) %>%
  mutate(X = ifelse(X=='(Intercept)', 'After', substring(X,17,23))) %>%
  left_join(counts, by=c('X'='Year')) %>%
  select(Source, Year = X, Count, p)

# split into some stuff pulling out the information and some stuff for plotting
# I think all we really need are quantiles

write.csv(pred, "data/attr_two_times.csv", row.names=FALSE)

library(forcats)
p <- pred %>% mutate(pC = p*Count, p = p*100) %>% gather("Version", "Value", p, pC) %>%
  mutate(Version = fct_recode(Version, Percent='p', Cases = 'pC'),
         Source = fct_relevel(Source, 'Poultry', 'Cattle', 'Sheep', 'Water', 'Other'),
         Year = fct_relevel(Year, 'Before'),
         Year = fct_recode(Year, `2005-2007`='Before', `2008-2016`='After'))

fig_width <- 960
fig_height <- 540

library(ggplot2)

#TODO: Flip this around

png("figures/05_attribution_both_b.png", width=fig_width, height=fig_height)
ggplot(p) +
  geom_violin(aes(Source, Value, fill=Source), scale='width') +
  facet_grid(Version~Year, scales = 'free_y') +
  theme_bw(base_size=20) +
  # theme(axis.text.x = element_text(hjust=c(0.2,0.5,0.5,0.5,0.8))) +
  xlab("") +
 # scale_x_discrete(limits=rev(levels(p$Source))) +
  scale_y_continuous(name = "Attributed human cases") +
#  coord_flip() +
  scale_fill_manual(values = c("plum4", "steelblue", "steelblue2", "brown", "green4"), guide=FALSE)

dev.off()

g1 <- ggplot(p %>% filter(Version == 'Percentage')) +
  geom_violin(aes(Source, Value, fill=Source), scale='width') +
  facet_grid(Year~.) +
  theme_bw(base_size=20) +
  # theme(axis.text.x = element_text(hjust=c(0.2,0.5,0.5,0.5,0.8))) +
  xlab("") +
  scale_x_discrete(limits=rev(levels(p$Source))) +
  scale_y_continuous(name = "Attributed human cases", labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values = c("plum4", "steelblue", "brown", "green4"), guide=FALSE)

g2 <- ggplot(p %>% filter(Version == 'Count')) +
  geom_violin(aes(Source, Value, fill=Source), scale='width') +
  facet_grid(Year~.) +
  theme_bw(base_size=20) +
  # theme(axis.text.x = element_text(hjust=c(0.2,0.5,0.5,0.5,0.8))) +
  xlab("") +
  scale_x_discrete(limits=rev(levels(p$Source))) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  scale_y_continuous(name = "Attributed human cases") +
  coord_flip() +
  scale_fill_manual(values = c("plum4", "steelblue", "brown", "green4"), guide=FALSE)

library(gridExtra)
grid.arrange(g1, g2, ncol=2)


