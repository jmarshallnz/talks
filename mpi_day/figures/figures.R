library(dplyr)

## STUFF FOR ALL PLOTS
alpha = function(col, alpha) { rgb(t(col2rgb(col)/255), alpha=alpha) }
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.03)
ax_col = "grey20"
fig_width = 10

###### NZ Cases of campy through time
cases <- read.csv("~/data/sa_report/data/dhb_cases.csv", stringsAsFactors = FALSE) %>%
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


###### ST dist on humans
library(forcats)
library(tidyr)
attr <- read.csv("~/data/sa_report/data/extract_attribution.csv")
sts <- attr %>% mutate(Source = fct_collapse(Source, Poultry = c("Supplier A", "Supplier B", "Supplier Other"),
                                            Water = "Environmental water"))

sts <- sts %>% group_by(ST, Source) %>% summarize(Count = n()) %>% spread(Source, Count, fill=0)
par(mar=c(3,3,0,1), mgp=c(2,.7,0), tck=-.015)
top20 = order(-sts$Human)[1:20]
barplot(c(sts$Human[top20], sum(sts$Human[-top20])), names=c(sts$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title(ylab="Cases", col.lab=ax_col)

###### ST dist on animals

sources <- c("Cattle", "Sheep", "Poultry", "Water")

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


###### ST dist by rurality

ur <- attr %>% filter(Source == "Human", !is.na(UR_bool)) %>% group_by(ST, UR_bool) %>% summarize(Count=n()) %>%
  spread(UR_bool, Count, fill=0)

par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,.7,0), tck=-.03)
top20 = match(sts$ST[order(-sts$Human)[1:20]], ur$ST)
cols = c(rep(alpha("steelblue",0.5), 20), "grey70")
barplot(c(ur$Urban[top20], sum(ur$Urban[-top20])), names=c(ur$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title("Urban", col.main=ax_col, cex.main=1.3)
title(ylab="Cases", col.lab=ax_col)
barplot(c(ur$Rural[top20], sum(ur$Rural[-top20])), names=c(ur$ST[top20], "Other"), col=cols, border=NA, yaxt="n", cex.names = 0.9, las=2, col.axis=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
title("Rural", col.main=ax_col, cex.main=1.3)
title(ylab="Cases", col.lab=ax_col)

###### Attribution results

#TODO: This should be urban/rural, but for some reason the fit seems very bad :(
attr_dyn = read.csv("~/data/sa_report/data/attribution_dynamic.csv", stringsAsFactors = FALSE)
source_map4 = read.csv("~/data/sa_report/source_maps/4_source.csv", stringsAsFactors = FALSE) %>% filter(Label != "Human") %>%
  dplyr::select(Number, Label) %>% unique %>% arrange(Number)

times = 1:max(attr_dyn$Month)
years = max(attr_dyn$Month)/12 - 1

urban = attr_dyn #[attr_dyn$UR_bool == "Rural",]
par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.03)
cols = c("plum4", "steelblue")
sources = c("Poultry", "Ruminants")
for (src in seq_along(sources)) {
  d = urban[urban$Source == sources[src],]
  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=d$Source[1], axes=FALSE, col.main=ax_col, ylab="", xlab="", xaxs="i", yaxs="i")
  polygon(c(times, rev(times)), c(d$ui[times], rev(d$li[times])), col=alpha(cols[src], 0.25), border=NA)
  polygon(c(times, rev(times)), c(d$yi[times], rev(d$xi[times])), col=alpha(cols[src], 0.35), border=NA)
  lines(times, d$mu[times], lwd=2, col=cols[src])
  axis(1, at=seq(0,144,by=12), labels=rep("", 13), col=ax_col, col.axis=ax_col, cex.axis=0.8)
  axis(2, at=seq(0,1,by=0.2), labels=paste0(seq(0,100,by=20),'%'), col=ax_col, col.axis=ax_col, cex.axis=0.8, las=1)
  mtext(2005:2016, side=1, col=ax_col, at=seq(0,132,by=12)+6, line=0.5)
}

####### Attribution results: Totals

# TODO: Ideally this would be UR results
attr_dyn = read.csv("~/data/sa_report/data/attribution_dynamic.csv", stringsAsFactors = FALSE)
source_map4 = read.csv("~/data/sa_report/source_maps/4_source.csv", stringsAsFactors = FALSE) %>% filter(Label != "Human") %>%
  dplyr::select(Number, Label) %>% unique %>% arrange(Number)

times = 1:max(attr_dyn$Month)
years = max(attr_dyn$Month)/12 - 1

urban = attr_dyn
par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.03)
cols = c("plum4", "steelblue")
sources = c("Poultry", "Ruminants")
for (src in seq_along(sources)) {
  d = urban[urban$Source == sources[src],]
  mu = d$mu * d$Count
  li = d$li * d$Count
  ui = d$ui * d$Count
  xi = d$xi * d$Count
  yi = d$yi * d$Count

  plot(NULL, xlim=range(times), ylim=c(0,40), type="n", main=d$Source[1], axes=FALSE, col.main=ax_col, ylab="", xlab="", xaxs="i", yaxs="i")
  polygon(c(times, rev(times)), c(ui[times], rev(li[times])), col=alpha(cols[src], 0.25), border=NA)
  polygon(c(times, rev(times)), c(yi[times], rev(xi[times])), col=alpha(cols[src], 0.35), border=NA)
  lines(times, mu[times], lwd=2, col=cols[src])
  axis(1, at=seq(0,144,by=12), labels=rep("", 13), col=ax_col, col.axis=ax_col, cex.axis=0.8)
  axis(2, at=seq(0,40,by=10), col=ax_col, col.axis=ax_col, cex.axis=0.8, las=1)
  title(ylab="Cases", col.lab=ax_col)
  mtext(2005:2016, side=1, col=ax_col, at=seq(0,132,by=12)+6, line=0.5)
}


####### Plot of Jing's work
cat <- read.table("~/data/sa_report/data/jing/cateI_capital.txt")
names(cat) <- c("Rurality", "Poultry", "Ruminants", "Water", "Other")

# drop down to quantiles and means
m <- cat %>% group_by(Rurality) %>% summarize_all(mean) %>% mutate(Variable='mu')
li <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.1)) %>% mutate(Variable='li')
ui <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.9)) %>% mutate(Variable='ui')
xi <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.3)) %>% mutate(Variable='xi')
yi <- cat %>% group_by(Rurality) %>% summarize_all(quantile, probs=c(0.7)) %>% mutate(Variable='yi')
d <- rbind(m, li, ui, xi, yi)
d <- d %>% gather('Source', 'Value', Poultry:Other) %>% spread(Variable, Value) %>% mutate_at(vars(mu:yi), function(x) { x*100 })

# TODO: Legend and colours same as above
library(ggplot2)
ggplot(d) + geom_ribbon(aes(x=Rurality, ymin=li, ymax=ui, fill=Source), alpha=0.25) + 
  geom_ribbon(aes(x=Rurality, ymin=xi, ymax=yi, fill=Source), alpha=0.35) + geom_line(aes(x=Rurality, y=mu, colour=Source)) +
  theme_bw() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), name="Percentage of cases")

####### Plot of ST-474 attribution, tidied up...

####### Plot of ST-474 tree/MDS ????
