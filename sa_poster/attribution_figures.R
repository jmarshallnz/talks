library(tidyr)
library(dplyr)
library(islandR)
library(ggplot2)
library(extrafont)

# attribution from Jing
d1 <- read.table("~/data/talks/sa_poster/data/jing/linearIsland_capitalF.txt")
names(d1) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")
d2 <- read.table("~/data/talks/sa_poster/data/jing/cateIsland_capitalF.txt")
names(d2) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")
d1$Model <- "Linear"
d2$Model <- "Categorical"
d1 <- rbind(d1, d2)
d1$Method <- "Island"

d3 <- read.table("~/data/talks/sa_poster/data/jing/linearDir_capitalF.txt")
names(d3) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")
d4 <- read.table("~/data/talks/sa_poster/data/jing/cateDir_capitalF.txt")
names(d4) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")
d3$Model <- "Linear"
d4$Model <- "Categorical"
d3 <- rbind(d3, d4)
d3$Method <- "Dirichlet"
d1 <- rbind(d1, d3)

order <- c("Poultry", "Ruminants", "Other", "Water")

# drop down to quantiles and means
m <- d1 %>% group_by(Model, Method, Rurality) %>% summarize_all(mean) %>% mutate(Variable='mu')
li <- d1 %>% group_by(Model, Method, Rurality) %>% summarize_all(quantile, probs=c(0.1)) %>% mutate(Variable='li')
ui <- d1 %>% group_by(Model, Method, Rurality) %>% summarize_all(quantile, probs=c(0.9)) %>% mutate(Variable='ui')
xi <- d1 %>% group_by(Model, Method, Rurality) %>% summarize_all(quantile, probs=c(0.3)) %>% mutate(Variable='xi')
yi <- d1 %>% group_by(Model, Method, Rurality) %>% summarize_all(quantile, probs=c(0.7)) %>% mutate(Variable='yi')
d <- rbind(m, li, ui, xi, yi)
d <- d %>% gather('Source', 'Value', Poultry:Other) %>% spread(Variable, Value) %>% mutate_at(vars(li:yi), function(x) { x*100 })
d$Source <- factor(d$Source, order)

#extrafont::loadfonts(device="pdf")
cairo_pdf("~/data/talks/sa_poster/figures/rural_attribution.pdf", width=15.5, height=11.5)
ggplot(d) + geom_ribbon(aes(x=Rurality, ymin=li, ymax=ui, fill=Source), alpha=0.3) +
  geom_ribbon(aes(x=Rurality, ymin=xi, ymax=yi, fill=Source), alpha=0.5) +
  geom_line(aes(x=Rurality, y=mu, colour=Source), size=1.5) +
  facet_grid(Model~Method) +
  theme_bw(base_size = 28) +
  scale_x_continuous(expand=c(0,0), name=NULL, labels=c('Highly Rural', rep('', 5), 'Highly Urban')) +
  scale_y_continuous(expand=c(0,0), name="Percentage of cases") +
  scale_color_manual(values = c("brown", "steelblue", "plum4", "green4")) +
  scale_fill_manual(values = c("brown", "steelblue", "plum4", "green4")) +
  theme(text = element_text(family="Cabin", colour="black"),
        legend.position = c(0.98,0.97),
        legend.justification = "right",
        legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(hjust=c(rep(-0.1,6),1.1), colour='black'),
        axis.text.y = element_text(colour='black'),
        panel.grid.major = element_line(colour='grey85'),
        panel.grid.minor = element_line(colour='grey85'),
        strip.background = element_rect(linetype="blank", fill="black"),
        strip.text.x = element_text(colour="grey80", face="bold", hjust=0, margin=margin(0.1,0.1,0.1,0.1,"inch")),
        strip.text.y = element_text(colour="grey80", face="bold", vjust=1, margin=margin(0.1,0.1,0.1,0.1,"inch"))) +
  guides(fill=guide_legend(title=NULL,nrow=1),
         color=guide_legend(title=NULL, nrow=1))
dev.off()






### NOW THE GENOTYPE DISTRIBUTION
library(tidyr)
library(dplyr)
library(islandR)
library(ggplot2)
library(ggjoy)

# testing shit
formula = Source ~ ST
non_primary = "Human"
data = manawatu

st = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            method = "dirichlet",
            data = manawatu)

set.seed(3)
st2 = st_fit(formula = Source ~ ST,
             non_primary = "Human",
             method="island",
             data = manawatu,
             sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
             iters=10000)
all.equal(st2, st2_old)

st2_old <- st2

# the island model will return more types in the sampling distribution than the dirichlet, so
# filter these out and reorder them accordingly (NOTE: This depends on knowing the structure
# of the sample_dist object - so may change in future)
st2$sampling_distribution <- st2$sampling_distribution[dimnames(st$sampling_distribution)[[1]],,]
st2$types <- st$types

# compute the sampling distributions by dividing by column sum
# (i.e. this is conditional on those we have observed so that the
# models are comparable)
ans <- lapply(1:100, function(x) { y = st$sampling_distribution[,,x]; z = data.frame(Iteration=x, ST=rownames(y), sweep(y, 2, colSums(y), '/'), stringsAsFactors = FALSE) })
ans <- do.call(rbind, ans)
ans2 <- lapply(1:100, function(x) { y = st2$sampling_distribution[,,x]; z = data.frame(Iteration=x, ST=rownames(y), sweep(y, 2, colSums(y), '/'), stringsAsFactors = FALSE) })
ans2 <- do.call(rbind, ans2)
ans$Model <- "Dirichlet"
ans2$Model <- "Island"
ans <- rbind(ans, ans2)
ans <- ans %>% gather("Source", "Probability", Other:Water)

# try breaking up the each Source by popular STs
count_source <- table(manawatu$Source, manawatu$ST)

top_st <- function(x, n) {
  names(x)[order(x, decreasing = TRUE)[seq_len(n)]]
}

top10 <- apply(count_source, 1, top_st, n=10)

get_subset <- function(source) {
  ans %>% filter(ST %in% top10[,source], Source == source) %>% mutate(ST = factor(ST, levels=rev(top10[,source])))
}

order <- c("Poultry", "Ruminants", "Other", "Water")
st_order <- paste(rep(order, each=10), as.character(top10[10:1,order]), sep="_")

plot_ans <- do.call(rbind, lapply(order, get_subset)) %>%
  mutate(ST_Source = factor(paste(Source, ST, sep="_"), st_order),
         Source = factor(Source, order))

# function for fixing the y axis labels
st_label <- function(x) sub("[^_]*_","",x )

cairo_pdf("~/data/talks/sa_poster/figures/genotype_dist.pdf", width=15.6, height=12)
ggplot(plot_ans) + geom_joy(aes(x=Probability, y=ST_Source, fill=Model), alpha=0.7, bandwidth=0.003, scale=1.2, size=0.1) +
  facet_wrap(~Source, scales="free_y") +
  scale_y_discrete(labels=st_label) +
  scale_x_continuous(limits=c(0,0.2), expand = c(0,0)) +
  scale_fill_manual(values = c("steelblue", "brown")) +
  theme_bw(base_size=28) +
  theme(text = element_text(family="Cabin", colour="black"),
        legend.position = c(0.907,0.085),
        legend.box.background = element_rect(),
        legend.margin = margin(6,12,12,12),
        legend.key.width = unit(30, "points"),
        legend.key.height = unit(30, "points"),
        plot.margin = unit(c(0.25,0.5,.25,.25), "inch"),
        axis.text.x = element_text(colour='black'),
        axis.text.y = element_text(colour='black'),
        panel.grid.major = element_line(colour='grey85'),
        panel.grid.minor = element_line(colour='grey85'),
        strip.background = element_rect(linetype="blank", fill="black"),
        strip.text.x = element_text(colour="grey80", face="bold", hjust=0, margin=margin(0.1,0.1,0.1,0.1,"inch")),
        strip.text.y = element_text(colour="grey80", face="bold", vjust=1, margin=margin(0.1,0.1,0.1,0.1,"inch"))) +
  guides(fill=guide_legend(title=NULL)) + ylab("Sequence Type")
dev.off()

