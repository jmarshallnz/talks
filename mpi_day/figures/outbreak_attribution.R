library(islandR)
library(dplyr)
library(readxl)
library(forcats)
library(pubmlst)
library(tidyr)

# read in data
manawatu_sources <- read.csv("data/sa_report/extract_attribution.csv")

# remap the sources
manawatu_sources <- manawatu_sources %>%
  filter(Source != "Human", Source != "Environmental water") %>%
#  filter(Source != "Environmental water", Source != "Wild_bird") %>%
#  mutate(SampledDate = as.Date(as.character(SampledDate))) %>%
#  filter(SampledDate >= "2014-01-01") %>%
  mutate(Source = as.character(fct_collapse(Source,
                               Poultry=c("Supplier A", "Supplier B",
                                         "Supplier Other", "Spent_hen",
                                         "Turkey", "Duck_poultry", "Otherpoultry"),
#                               Ruminants=c("Sheep", "Cattle"),
                               Other=c("Cat_dog_pet", "Pig", "Water_bird_wild", "Wild_bird_other"))))

# HUMAN DATA FROM PAHIATUA
pahiatua <- data.frame(ST=c(61,190,190,190,3793,50,190,50,61))
cases <- pahiatua %>% left_join(pubmlst) %>% select(-CC, -Coli)
cases$Source <- "Human"

dat <- rbind(cases, manawatu_sources %>% select(one_of(names(cases))))
dat$Source <- factor(dat$Source)
dat <- data.frame(dat)

# estimate sampling distribution on sources
st = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = dat,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC)

# attribute human cases
mod = attribution(ST ~ 1, st, data=subset(dat, Source == "Human"))

summary(mod)
quantiles <- predict(mod, newdata=NULL, FUN=quantile, c(0.025, 0.5, 0.975))
attribution <- quantiles %>%
  select(-X) %>% mutate(p = round(p*100, 1)) %>%
  spread(quantile, p)

restrict <- quantiles %>% spread(quantile, p) %>% rename(lc = `2.5%`, uc=`97.5%`)
posterior = predict(mod, FUN=identity)

filtered <- posterior %>% left_join(restrict) %>%
  filter(p > lc, p < uc)

write.csv(filtered, "data/pahiatua.csv", row.names=FALSE)

### RAW MILK OUTBREAK
rawmilk <- read.csv("data/outbreak_rawmilk.csv") %>% select(ST) %>%
  left_join(pubmlst) %>% select(-CC, -Coli) %>% mutate(Source = "Human")

dat <- rbind(rawmilk, manawatu_sources %>% select(one_of(names(cases))))
dat$Source <- factor(dat$Source)
dat <- data.frame(dat)


# estimate sampling distribution on sources
st = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = dat,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC)

# attribute human cases
mod = attribution(ST ~ 1, st, data=subset(dat, Source == "Human"))

summary(mod)
quantiles <- predict(mod, newdata=NULL, FUN=quantile, c(0.025, 0.5, 0.975))
attribution <- quantiles %>%
  select(-X) %>% mutate(p = round(p*100, 1)) %>%
  spread(quantile, p)

restrict <- quantiles %>% spread(quantile, p) %>% rename(lc = `2.5%`, uc=`97.5%`)
posterior = predict(mod, FUN=identity)

filtered <- posterior %>% left_join(restrict) %>%
  filter(p > lc, p < uc)

write.csv(filtered, "data/rawmilk.csv", row.names=FALSE)


# plot the attribution
outbreaks <- rbind(read.csv("data/pahiatua.csv") %>% mutate(Outbreak="Pahiatua"),
                   read.csv("data/rawmilk.csv") %>% mutate(Outbreak="Raw milk")) %>%
  mutate(Source = fct_relevel(Source, "Poultry", "Cattle", "Sheep", "Other"))

png("figures/15_outbreak_attribution.png", width=fig_width, height=fig_height)
ggplot(outbreaks) +
  geom_violin(aes(Source, p, fill=Source), scale='width') +
  theme_bw(base_size=20) +
  facet_grid(. ~ Outbreak) +
  scale_fill_manual(values = c("plum4", "steelblue", "steelblue2", "brown"), guide = "none") +
  xlab("") +
  coord_flip() +
  scale_y_continuous(name = "Attribution of human cases", labels = scales::percent) +
  scale_x_discrete(limits=rev(levels(outbreaks$Source)))
dev.off()
