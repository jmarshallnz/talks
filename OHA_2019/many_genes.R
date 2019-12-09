
base_path <- "~/data/OneDrive/projects/spencer_attribution/sacnz/output/"

attr1 <- read_csv(file.path(base_path, "per_source/attribution.csv")) %>%
  mutate(Source = paste0("Source", as.numeric(as.factor(Source))), Model="PerSource")
attr2 <- read_csv(file.path(base_path, "common/attribution.csv")) %>%
  mutate(Source = paste0("Source", as.numeric(as.factor(Source))), Model="Common")

best_attr <- bind_rows(attr1, attr2) %>% filter(str_detect(Genes, "Best"), str_detect(Data, "Full")) %>%
  extract(Genes, into='Genes', regex="([0-9]+)", convert=TRUE) %>%
  spread(quantile, value) %>% select(Model, Source, Genes, p=`0.5`) %>%
  nest(data = c(Genes, p)) %>% mutate(model = map(data, ~loess(p ~ Genes, data=.))) %>%
  mutate(fits = map2(model, data, predict)) %>%
  unnest(c(data, fits)) %>% select(-model)

write_csv(best_attr, "data/many_genes.csv")
ggplot(best_attr, aes(x=Genes, y=fits, fill=Source, col=Source)) +
  geom_ribbon(aes(ymin=fits-0.03, ymax = fits+0.03), alpha=0.4, col='transparent') +
  geom_line(lwd=1) +
  scale_y_continuous("Attribution of human cases", labels = scales::percent, expand=c(0,0), limits = c(-1,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_x_continuous(limits = c(5,100), expand=c(0,0)) +
  guides(col = 'none', fill='none') + facet_wrap(~Model)

attr <- read_csv(file.path(base_path, "common/attribution.csv")) %>%
  mutate(Source = paste0("Source", as.numeric(as.factor(Source))))
best_attr <- attr %>% filter(str_detect(Genes, "Best"), str_detect(Data, "Full")) %>%
  extract(Genes, into='Genes', regex="([0-9]+)", convert=TRUE) %>%
  spread(quantile, value) %>% select(Source, Genes, p=`0.5`) %>%
  nest(data = c(Genes, p)) %>% mutate(model = map(data, ~loess(p ~ Genes, data=.))) %>%
  mutate(fits = map2(model, data, predict)) %>%
  unnest(c(data, fits))

ggplot(best_attr, aes(x=Genes, y=p, fill=Source, col=Source)) +
  geom_smooth(level=0.999999) +
  scale_y_continuous("Attribution of human cases", labels = scales::percent, expand=c(0,0), limits = c(-1,1)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_x_continuous(limits = c(5,100), expand=c(0,0)) +
  guides(col = 'none', fill='none')
