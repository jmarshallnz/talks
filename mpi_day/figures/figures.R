library(dplyr)
library(forcats)
library(tidyr)
library(ggplot2)
library(mgcv)

data_folder <- "data/sa_report"
fig_width <- 960
fig_height <- 540

## STUFF FOR ALL PLOTS
alpha = function(col, alpha) { rgb(t(col2rgb(col)/255), alpha=alpha) }
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.03)
ax_col = "grey20"

###### NZ Cases of campy through time
cases <- read.csv(file.path(data_folder, "dhb_cases.csv"), stringsAsFactors = FALSE) %>%
  mutate(Population = PopulationInterpolated)

# filter out where we don't have a population
cases <- cases %>% filter(!is.na(Population))

nz_cases <- cases %>% select(Year, Month, Count, Population) %>%
  group_by(Year, Month) %>% na.omit %>%
  summarise(Count=sum(Count), Population=sum(Population), Rate=Count/Population*100000*12) %>%
  mutate(MonthNo = match(Month, month.name),
         Date = sprintf('%04d-%02d-%02d', Year, MonthNo, 1),
         Month = month.abb[MonthNo]) %>%
  arrange(Date)

# do the plot
nz_cases$Date <- as.Date(nz_cases$Date)
nz_cases$Month <- as.factor(nz_cases$Month)

mod.gam = gam(Rate ~ s(as.numeric(Date)) + Month, data=nz_cases)
co = coef(mod.gam)
avg_month = mean(co[grepl("Month",names(co))])

png("figures/01_nzcases.png", width = fig_width, height = fig_height)
par(mar=c(3,3.5,0.5,1), mgp=c(2,.7,0), tck=-.015)
plot(Rate ~ Date, data=nz_cases, type="l", col=ax_col, axes=FALSE, xaxs="i", yaxs="i", ylim=c(0,600), xlab="", ylab="", xlim=as.Date(c("2006-01-01", "2017-01-01")))
y = predict(mod.gam, data.frame(Date=as.numeric(nz_cases$Date), Month="Apr"), se.fit=TRUE)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2005:2017,"-01-01")), labels=rep("",13))
mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
y_fit = y$fit + avg_month
polygon(c(nz_cases$Date, rev(nz_cases$Date)), c(y_fit+y$se.fit, rev(y_fit - y$se.fit)), col=alpha("steelblue", 0.3), border=NA)
lines(nz_cases$Date, y_fit, col="steelblue", lwd=2)
title(ylab="Cases per 100,000", col.lab=ax_col, line=2.5)
dev.off()

###### ST dist on humans
attr <- read.csv(file.path(data_folder, "extract_attribution.csv"))
sts <- attr %>% mutate(Source = fct_collapse(Source, Poultry = c("Supplier A", "Supplier B", "Supplier Other"),
                                            Water = "Environmental water"))

sts <- sts %>% group_by(ST, Source) %>% summarize(Count = n()) %>% spread(Source, Count, fill=0)

png("figures/02_st_human.png", width = fig_width, height = fig_height)
par(mar=c(3,3,0,1), mgp=c(2,.7,0), tck=-.015)
top20 = order(-sts$Human)[1:20]
barplot(c(sts$Human[top20], sum(sts$Human[-top20])), names=c(sts$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title(ylab="Cases", col.lab=ax_col)
dev.off()

###### ST dist on animals

sources <- c("Cattle", "Sheep", "Poultry", "Water")

png("figures/03_st_sources.png", width = fig_width, height = fig_height)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.02)
for (s in sources) {
  cols = c(rep(alpha("steelblue",0.5), 20), "grey70")
  barplot(c(sts$Human[top20])/sum(sts$Human)*100, names=sts$ST[top20], col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
  cols = c(rep(alpha("plum4",0.8), 21))
  scale <- sum(sts$Human)
  barplot((sts[top20, s] %>% unlist) / sum(sts[,s])*100, names=rep("", 20), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col, add=TRUE)
  #  x = unlist(par()["mfg"])[2]
  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  #  text(ifelse(x == 1, 5, 0.2), 70, st, col = ax_col, cex = 2, adj=c(ifelse(x == 1, 1, 0), 0.5))
  title(ylab="Percentage", col.lab=ax_col)
  title(main = s)
}
dev.off()


###### ST dist by rurality

ur <- attr %>% filter(Source == "Human", !is.na(UR_bool)) %>% group_by(ST, UR_bool) %>% summarize(Count=n()) %>%
  spread(UR_bool, Count, fill=0)

png("figures/04_st_urban_rural.png", width = fig_width, height = fig_height)
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,.7,0), tck=-.03)
top20 = match(sts$ST[order(-sts$Human)[1:20]], ur$ST)
cols = c(rep(alpha("steelblue",0.7), 20), "grey50")
barplot(c(ur$Urban[top20], sum(ur$Urban[-top20])), names=c(ur$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title("Urban", col.main=ax_col, cex.main=1.3)
title(ylab="Cases", col.lab=ax_col)
barplot(c(ur$Rural[top20], sum(ur$Rural[-top20])), names=c(ur$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title("Rural", col.main=ax_col, cex.main=1.3)
title(ylab="Cases", col.lab=ax_col)
dev.off()

###### Attribution results

#TODO: This should be urban/rural, but for some reason the fit seems very bad :(
attr_dyn = read.csv(file.path(data_folder, "attribution_dynamic_urban_rural_TEST_PRIORS.csv"), stringsAsFactors = FALSE)
source_map4 = read.csv(file.path(data_folder, "4_source.csv"), stringsAsFactors = FALSE) %>% filter(Label != "Human") %>%
  dplyr::select(Number, Label) %>% unique %>% arrange(Number)

times = 1:max(attr_dyn$Month)
years = max(attr_dyn$Month)/12 - 1
attr_dyn <- attr_dyn %>% mutate(Year = 2005 + (Month-1) %/% 12,
                    Month = (Month-1) %% 12 + 1,
                    Date = as.Date(sprintf('%04d-%02d-01', Year, Month)),
                    UR_bool = fct_relevel(UR_bool, "Urban", "Rural"))

png("figures/05_attribution_proportion.png", width = fig_width, height = fig_height)
ggplot(attr_dyn %>% filter(Source %in% c('Poultry', 'Ruminants'))) + 
  geom_ribbon(aes(x=Date, ymin=li, ymax=ui, fill=Source), alpha=0.4) +
  geom_ribbon(aes(x=Date, ymin=xi, ymax=yi, fill=Source), alpha=0.7) +
  geom_line(aes(x=Date, y=mu)) +
  facet_grid(Source ~ UR_bool) +
  theme_bw(base_size=20) +
  scale_fill_manual(values = c("plum4", "steelblue"), guide = FALSE) +
  scale_x_date(expand=c(0, 0)) +
  scale_y_continuous(name = "Attributed human cases", labels = scales::percent)
dev.off()

#### TOTALS

png("figures/06_attribution_totals.png", width = fig_width, height = fig_height)
ggplot(attr_dyn %>% filter(Source %in% c('Poultry', 'Ruminants'))) + 
  geom_ribbon(aes(x=Date, ymin=li*Count, ymax=ui*Count, fill=Source), alpha=0.4) +
  geom_ribbon(aes(x=Date, ymin=xi*Count, ymax=yi*Count, fill=Source), alpha=0.7) +
  geom_line(aes(x=Date, y=mu*Count)) +
  facet_grid(Source ~ UR_bool) +
  theme_bw(base_size=20) +
  scale_fill_manual(values = c("plum4", "steelblue"), guide = FALSE) +
  scale_x_date(expand=c(0, 0)) +
  scale_y_continuous(name="Attributed human cases")
dev.off()


####### Plot of Jing's work
if (0) {
  cat <- read.table(file.path(data_folder, "cateI_capital.txt"))
  names(cat) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")
  
  # drop down to quantiles and means
  m <- cat %>% group_by(Rurality) %>% summarize_all(mean) %>% mutate(Variable='mu')
  li <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.1)) %>% mutate(Variable='li')
  ui <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.9)) %>% mutate(Variable='ui')
  xi <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.3)) %>% mutate(Variable='xi')
  yi <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.7)) %>% mutate(Variable='yi')
  d <- rbind(m, li, ui, xi, yi)
  d <- d %>% gather('Source', 'Value', Poultry:Other) %>% spread(Variable, Value) %>% mutate_at(vars(li:yi), function(x) { x*100 })
  
  write.csv(d, file.path(data_folder, 'ur_categorical.csv'), row.names = FALSE)
}

d <- read.csv(file.path(data_folder, 'ur_categorical.csv'))
d$Source <- factor(d$Source, c('Poultry', 'Ruminants', 'Water', 'Other'))

prior_labeller <- function(x) {
  x <- ifelse(x == "0", "No expert reliance", ifelse(x == "2", "Strong expert reliance", ""))
}

png("figures/07_attribution_rurality.png", width = fig_width, height = fig_height)
ggplot(d) + geom_ribbon(aes(x=Rurality, ymin=li, ymax=ui, fill=Source), alpha=0.25) + 
  geom_ribbon(aes(x=Rurality, ymin=xi, ymax=yi, fill=Source), alpha=0.35) + geom_line(aes(x=Rurality, y=mu, colour=Source)) +
  theme_bw(base_size = 20) +
  scale_x_continuous(expand=c(0,0), labels=c('Highly Rural', rep('', 5), 'Highly Urban'), name='') +
  scale_y_continuous(expand=c(0,0), name="Percentage of cases") +
  scale_color_manual(values = c("plum4", "steelblue", "brown", "green4")) +
  scale_fill_manual(values = c("plum4", "steelblue", "brown", "green4")) +
  theme(legend.position = c(0.95,0.95),
        legend.justification = "right",
        legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(hjust=c(rep(-0.1,6),1.1))) +
  guides(fill=guide_legend(title=NULL,nrow=1),
         color=guide_legend(title=NULL, nrow=1))
dev.off()
