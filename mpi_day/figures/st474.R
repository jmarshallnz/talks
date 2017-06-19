library(readxl)
library(islandR)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(igraph)
library(cluster)

source('figures/tree_stuff.R')

####### Plotting diversity of ST474

wgmlst <- read.table('data/474/allele_profiles.txt')
names(wgmlst)[1] <- "genome"

# remove columns that are useless
v_col <- apply(wgmlst[,-1],2,var)

wgmlst <- wgmlst[, c(1,which(v_col > 0)+1)]

# add an ST column, being the pasted copy of everything else
wgmlst$ST <- as.numeric(as.factor(apply(wgmlst[,-1],1,paste,collapse='_')))

iso <- read_excel('data/474/st474_530isolates.xlsx') %>%
  mutate(Date = as.Date(ifelse(is.na(SampledDate), as.character(TestedDate), as.character(SampledDate)))) %>%
  select(genome, source=SA_model_source, Date) %>%
  mutate(source = ifelse(source == 'Supplier A', ifelse(Date < '2008-01-01', 'SuppA_before', 'SuppA_after'), source))

wgmlst <- wgmlst %>% left_join(iso) %>%
  mutate(source = fct_collapse(source,
                               `Poultry A 2005-2007` = 'SuppA_before',
                               `Poultry A 2008-2015` = 'SuppA_after',
                               `Poultry Other` = c("Supplier_other", "Supplier B", "Spent_hen"),
                               Ruminant = c("Cattle", "Sheep"),
                               Other = c("Cat_dog_pet", "Pig", "Wild_bird_other"),
                               Water = "Environmental water",
                               Human = "Human")) %>%
  mutate(source = fct_relevel(source, 'Poultry A 2005-2007', 'Poultry A 2008-2015', 'Poultry Other', 'Ruminant', 'Other', 'Water', 'Human'))

# filter out stuff we don't want
final <- wgmlst %>% filter(!is.na(source)) # %>% filter(source %in% c('A_after', 'A_before', 'Environmental water', 'Human', 'Supplier_other', 'Cattle', 'Sheep'))

# TODO: Try and do some plotting of this shit
loci <- final %>% select(V2:V1544) %>% mutate_all(factor)

d <- daisy(loci)

o <- MASS::isoMDS(d, maxit=100)
plot(o$points)

cols <- c("pink", "plum4", "purple", "steelblue", "brown", "green4", "orange")
g <- graph_from_adjacency_matrix(as.matrix(d), mode='undirected', weighted=TRUE)
t <- mst(g)
set.seed(7)
layout = layout_nicely(t)
layout2 <- norm_coords(layout, xmin = -1.4, xmax=0.6)

# try reducing the tree to remove all the ones we don't want
wch <- which(final$source != "Human")
t2 <- reduce_tree(t, wch)
layout3 = layout2[match(V(t2)$name, V(t)$name),]

png("figures/08_474_tree_sources.png", width=fig_width, height=fig_height)
par(mar=c(0,0,0,0))
plot(t2, vertex.label=NA, vertex.size=4, vertex.color = cols[final$source[wch]], layout = layout3, rescale=FALSE)
legend('topright', legend=levels(final$source)[1:6], col='black', pch=21, pt.bg = cols, bty='n', cex=1.7)
dev.off()

png("figures/09_474_tree_humans.png", width=fig_width, height=fig_height)
par(mar=c(0,0,0,0))
vsize <- rep(3,nrow(final))
vsize[wch] <- 4
plot(t, vertex.label=NA, vertex.size=vsize, vertex.color = cols[final$source], layout = layout2, rescale = FALSE)
legend('topright', legend=levels(final$source), col='black', pch=21, pt.bg = cols, bty='n', cex=1.7)
dev.off()

# ATTRIBUTION

r <- read.csv("data/474/attribution.csv") %>%
  mutate(Source = fct_relevel(Source, 'Poultry A 2005-2007', 'Poultry A 2008-2015', 'Poultry Other', 'Ruminant', 'Other', 'Water'))

png("figures/10_474_attribution.png", width=fig_width, height=fig_height)
r <- r %>% mutate(Col = cols[Source])
ggplot(r) +
  geom_violin(aes(Source, p, fill=Source), scale='width') +
  facet_grid(.~X) +
  theme_bw(base_size=20) +
  theme(panel.spacing = unit(0.8, "cm"), plot.margin=unit(c(0.5,1,0.5,0.5), 'cm'),
        axis.text.x = element_text(hjust=c(rep(-0.1,6),1.1))) +
  xlab("") +
  scale_x_discrete(limits=rev(levels(r$Source))) +
  scale_y_continuous(name = "Attributed human cases", labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values = cols, guide=FALSE)
dev.off()

