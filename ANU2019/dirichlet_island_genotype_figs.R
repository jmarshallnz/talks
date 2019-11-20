library(tidyverse)
library(ggridges)

load("~/data/R/dirichlet_island/fig_genotype_dist.Rdata")

# which STs to plot


final %>% gather(Model, P, Dirichlet, Island) %>%
  group_by(ST, Source, Model) %>% summarise(m = median(P), lci = quantile(P, 0.025),
                                                    uci = quantile(P, 0.975), ucl = quantile(P, 0.8),
                                                    lcl = quantile(P, 0.2)) %>% ungroup() -> foo

write_csv(foo, "data/genotypes.csv")

foo <- read_csv("data/genotypes.csv") %>% mutate(Source = factor(Source, levels = c("Poultry", "Ruminants", "Other", "Water")))

common_sts <- c(474, 45, 50, 53, 48, 61, 190, 42, 354, 520)
plot_common <- foo %>% filter(ST %in% common_sts) %>%
  mutate(ST = factor(ST, levels=common_sts, labels = paste0("ST-", common_sts)))

sts <- c(403, 2343, 2026, 474)
plot_dat <- foo %>% filter(ST %in% sts) %>% 
  mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)))

ggplot(plot_dat) + geom_linerange(aes(x=Source, ymin = lci, ymax=uci, group=Model), position=position_dodge(0.5)) + 
  geom_linerange(aes(x=Source, ymin = lcl, ymax = ucl, col=Model), size=2, position=position_dodge(0.5)) +
  geom_point(aes(x=Source, y=m, fill=Model), position=position_dodge(0.5), shape=21, size=2) +
  facet_wrap(~ST, nrow=2) +
  xlab("") +
  scale_y_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("steelblue2", "brown")) +
  scale_colour_manual(values = c("steelblue2", "brown")) +
  theme_bw() +
  theme(legend.position = c(0.92,0.895),
        legend.title = element_blank())

ggplot(plot_common) +  geom_boxplot(aes(x=Source, y=P, fill=Model), alpha=0.7, coef=5, size=0.1) +
  facet_wrap(~ST, nrow=2) +
  xlab("") +
  scale_y_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("steelblue2", "brown")) +
  theme_bw() +
  theme(legend.position = c(0.92,0.895),
        legend.box.background = element_rect(),
        legend.margin = margin(1.5,3,3,3),
        legend.title = element_blank())
#        axis.text.x = element_text(hjust=c(-0.1,rep(0.5, 3), 1.1)))


plot_dat <- final %>% filter(ST %in% sts) %>%
  gather(Model, P, Dirichlet, Island) %>% 
  filter(Model == "Island") %>% ungroup %>% mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)))

#pdf("fig_genotype_dist.pdf", width=7, height=4)
ggplot(plot_dat) + geom_density_ridges(aes(x=P, y=Source), fill='darkred', alpha=0.7, size=0.1, bandwidth=0.01) +
  facet_wrap(~ST, ncol=2) +
  ylab("") +
  scale_x_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  #  scale_fill_manual(values = c("steelblue2", "brown")) +
  theme_bw() +
  theme(legend.position = c(0.92,0.895),
        legend.box.background = element_rect(),
        legend.margin = margin(1.5,3,3,3),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust=c(-0.1,rep(0.5, 3), 1.1)))

ggplot(plot_dat) + geom_violin(aes(x=Source, y=P), fill='darkred', alpha=0.7) +
  facet_wrap(~ST, ncol=2) +
  ylab("") +
  scale_y_continuous(name = "P(Source | ST)", limits=c(0,1), expand = c(0,0)) +
  #  scale_fill_manual(values = c("steelblue2", "brown")) +
  theme_bw() +
  theme(legend.position = c(0.92,0.895),
        legend.box.background = element_rect(),
        legend.margin = margin(1.5,3,3,3),
        legend.title = element_blank(),
        axis.text.x = element_text(hjust=c(-0.1,rep(0.5, 3), 1.1)))

dev.off()
