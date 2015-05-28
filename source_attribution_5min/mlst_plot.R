# read the data in and do a hierarchical clustering tree thingee
# (probably using average linkage) across all types, and then
# histograms of frequency per type, ordered by the order in the
# clustering.

# then show the clustering with the histograms next to them.
# We can perhaps consider doing an animated or multi-slide one
# for pre/post intervention on the humans.

# or maybe seasonality?

library(dplyr)
library(cluster)
library(reshape2)

dat <- read.csv("E:/campy_db/reports/attribution_isolate_data/20150329.csv", stringsAsFactors=F)

# combine some sources together
dat <- dat %>% mutate(Source = ifelse(substring(Source,1,8) == "Supplier", "Poultry", Source))
dat <- dat %>% mutate(Source = ifelse(Source == "Environmental water", "Water", Source))
dat <- dat %>% mutate(Source = ifelse(Source != "Cattle" & Source != "Sheep" & Source != "Poultry" & Source != "Human" & Source != "Water", "Other", Source))

# filter out any that we've only ever observed once, as they're a bit boring
tab <- table(dat$ST)
rownames(tab)[which(tab > 1)]
frequent_sts <- as.numeric(names(which(table(dat$ST) > 1)))

dat <- dat %>% filter(ST %in% frequent_sts)

mlst_cols <- c("ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")

mlst_types <- dat %>% select(one_of(mlst_cols), ST) %>% unique %>% arrange(ST)
rownames(mlst_types) <- mlst_types$ST
mlst_types$ST <- NULL

# compute distance matrix and hierarchical clustering
for (i in 1:ncol(mlst_types))
  mlst_types[,i] <- as.factor(mlst_types[,i])
mlst_dist <- daisy(mlst_types)

mlst_hc <- hclust(mlst_dist, method="average")

plot(mlst_hc, hang=-1)

# compute barplots of frequencies per source
totals <- dat %>% group_by(Source, ST) %>% summarize(Total=n())
totals <- dcast(totals, ST ~ Source)
totals[is.na(totals)] <- 0

totals <- totals[mlst_hc$order,]
rownames(totals) <- totals$ST
totals$ST <- NULL
totals <- sweep(totals, 2, colSums(totals), FUN="/")
totals

# rearrange
totals <- totals %>% select(Poultry, Cattle, Sheep, Water, Other, Human)
totals <- totals / 0.252 #max(totals)


# function for modifying colours to make them transparent
alpha <- function(col, a)
{
  rgb(t(col2rgb(col)/255), alpha=a)
}


# make the barplots drip down onto a colour bar thingee
vals <- sweep(totals, 2, apply(totals, 2, max), FUN="/")

o <- c(1,2,3,4,6,8,8)
col = c("#FF7F00","#CF0000","#004FCF", "#009F9F","#8F006F","#9F5F3F","#FFAFAF", "#000000")
col <- col[o]

# read in the attribution data and draw a diagram on the right
#attribution <- data.frame(mean=c(0.5,0.3,0.2,0.2), lci=c(0.3,0.1,0.1,0.1), uci=c(0.7,0.4,0.3,0.3))
attribution <- read.csv("sa_presentation_out.csv")

all <- attribution %>% filter(name == "All")
png("figures/mlst1.png", width=720, height=440)
par(mai=c(3,0.5,0,0.5))
plot(mlst_hc, hang=-1, ylim=c(-1,1), labels=FALSE, xaxt="n", xlab="", sub="", main="", axes=FALSE, yaxt="n", ylab="", ann=FALSE)
for (i in 1:6) {
  barplot(rep(0.17, nrow(totals)), add=TRUE, offset=-0.18*i, col=alpha(col[i], vals[,i]*0.9+0.1), border=NA, space=0, axes=FALSE)
  barplot(-totals[,i]*0.12, add=TRUE, offset=-0.18*i+0.17, col=col[i], border=NA, space=0, axes=FALSE)
  text(-1, -0.18*i+0.09, names(totals)[i],adj=c(1,0.5), xpd=NA, cex=1.2)
}
for (i in seq_along(all$source)) {
  x <- nrow(mlst_types)+1
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$uci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$mean[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$lci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
}
dev.off()

# let's do another one without other + water
png("figures/mlst2.png", width=720, height=400)
par(mai=c(2,0.5,0.25,0.5))
plot(mlst_hc, hang=-1, ylim=c(-1,1), labels=FALSE, xaxt="n", xlab="", sub="", main="", axes=FALSE, yaxt="n", ylab="", ann=FALSE)
wch <- c(1:3,6)
for (i in seq_along(wch)) {
  barplot(rep(0.17, nrow(totals)), add=TRUE, offset=-0.18*i, col=alpha(col[wch[i]], vals[,wch[i]]*0.9+0.1), border=NA, space=0, axes=FALSE)
  barplot(-totals[,wch[i]]*0.12, add=TRUE, offset=-0.18*i+0.17, col=col[wch[i]], border=NA, space=0, axes=FALSE)
  text(-1, -0.18*i+0.09, names(totals)[wch[i]],adj=c(1,0.5), xpd=NA, cex=1.2)
}
for (i in 1:3) {
  x <- nrow(mlst_types)+1
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$uci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$mean[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  rect(x, -0.18*i, nrow(mlst_types)+1 + 15*all$lci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
}
#x <- order(-rowSums(totals[,wch]))[1:10]
#labs <- data.frame(x = x, lab = rownames(totals)[x], v = rowSums(totals[x,wch])) %>% arrange(x)
#labs_odd <- labs[seq(1,nrow(labs),by=2),]
#labs_even <- labs[seq(2,nrow(labs),by=2),]
#text(labs_odd$x-0.5, -0.18*4-0.01, labels=labs_odd$lab, cex=0.6, xpd=NA, adj=c(0.5,1))
#text(labs_even$x-0.5, -0.18*4-0.03, labels=labs_even$lab, cex=0.6, xpd=NA, adj=c(0.5,1))
dev.off()

# now figure out contribution pre- and post- intervention etc.

totals <- dat %>% group_by(Source, ST) %>% summarize(Total=n(), Urban=sum(UR_bool == "Urban", na.rm=T),
                                                     Rural=sum(UR_bool == "Rural", na.rm=T),
                                                     Before=sum(Intervention == "before", na.rm=T),
                                                     After=sum(Intervention == "after", na.rm=T),
                                                     UrbanBefore=sum(UR_bool == "Urban" & Intervention == "before", na.rm=T),
                                                     UrbanAfter=sum(UR_bool == "Urban" & Intervention == "after", na.rm=T),
                                                     RuralBefore=sum(UR_bool == "Rural" & Intervention == "before", na.rm=T),
                                                     RuralAfter=sum(UR_bool == "Rural" & Intervention == "after", na.rm=T))

vars <- names(totals)[4:ncol(totals)]
for (var in vars) {
  attr <- attribution %>% filter(name == var)

  total.Select <- dcast(totals, ST ~ Source, value.var=var) %>% select(ST, Human)
  total.Sources <- dcast(totals, ST ~ Source, value.var="Total") %>% select(Poultry, Cattle, Sheep, Water, Other, AllHuman=Human)
  total.All <- cbind(total.Sources, total.Select)
  total.All[is.na(total.All)] <- 0
  
  # reorder
  rownames(total.All) <- total.All$ST
  total.All$ST <- NULL
  total.All <- total.All[mlst_hc$order,]
  
  # scale
  total.All <- sweep(total.All, 2, colSums(total.All), FUN="/")
#  cat(max(total.All), "\n")
  total.All <- total.All / 0.252 #max(total.All)
  
  vals <- sweep(total.All, 2, apply(total.All, 2, max), FUN="/")
  
  png(paste0("figures/mlst_", var, ".png"), width=720, height=440)
  par(mai=c(3,0.5,0,0.5))
  plot(mlst_hc, hang=-1, ylim=c(-1,1), labels=FALSE, xaxt="n", xlab="", sub="", main="", axes=FALSE, yaxt="n", ylab="", ann=FALSE)
  wch <- c(1:5,7)
  for (i in seq_along(wch)) {
    barplot(rep(0.17, nrow(total.All)), add=TRUE, offset=-0.18*i, col=alpha(col[wch[i]], vals[,wch[i]]*0.9+0.1), border=NA, space=0, axes=FALSE)
    barplot(-total.All[,wch[i]]*0.12, add=TRUE, offset=-0.18*i+0.17, col=col[wch[i]], border=NA, space=0, axes=FALSE)
    text(-1, -0.18*i+0.09, names(total.All)[wch[i]],adj=c(1,0.5), xpd=NA, cex=1.2)
  }
  for (i in 1:5) {
    x <- nrow(mlst_types)+1
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$uci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$mean[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$lci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  }
  dev.off()

}



poultry <- read.csv("sa_poultry_presentation_out.csv")

vars <- c("Before", "After")
for (var in vars) {
  attr <- poultry %>% filter(name == var)

  total.Select <- dcast(totals, ST ~ Source, value.var=var) %>% select(ST, Poultry, Human)
  total.Sources <- dcast(totals, ST ~ Source, value.var="Total") %>% select(AllPoultry=Poultry, Cattle, Sheep, AllHuman=Human, Other, Water)
  total.All <- cbind(total.Sources, total.Select)
  total.All[is.na(total.All)] <- 0
  total.All <- total.All %>% select(ST, Poultry, Cattle, Sheep, Water, Other, Human, AllPoultry, AllHuman)
  
  # reorder
  rownames(total.All) <- total.All$ST
  total.All$ST <- NULL
  total.All <- total.All[mlst_hc$order,]
  
  # scale
  total.All <- sweep(total.All, 2, colSums(total.All), FUN="/")
  total.All <- total.All / 0.252 #max(total.All)
  
  vals <- sweep(total.All, 2, apply(total.All, 2, max), FUN="/")
  
  png(paste0("figures/mlst_poultry_", var, ".png"), width=720, height=440)
  par(mai=c(3,0.5,0,0.5))
  plot(mlst_hc, hang=-1, ylim=c(-1,1), labels=FALSE, xaxt="n", xlab="", sub="", main="", axes=FALSE, yaxt="n", ylab="", ann=FALSE)
  wch <- c(1:6)
  for (i in seq_along(wch)) {
    barplot(rep(0.17, nrow(total.All)), add=TRUE, offset=-0.18*i, col=alpha(col[wch[i]], vals[,wch[i]]*0.9+0.1), border=NA, space=0, axes=FALSE)
    barplot(-total.All[,wch[i]]*0.12, add=TRUE, offset=-0.18*i+0.17, col=col[wch[i]], border=NA, space=0, axes=FALSE)
    text(-1, -0.18*i+0.09, names(total.All)[wch[i]],adj=c(1,0.5), xpd=NA, cex=1.2)
  }
  for (i in 1:5) {
    x <- nrow(mlst_types)+1
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$uci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$mean[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
    rect(x, -0.18*i, nrow(mlst_types)+1 + 15*attr$lci[i], -0.18*i + 0.17, col=alpha(col[i],0.3), border=NA, xpd=NA)
  }
#  text(1:nrow(mlst_types), 0.18*4, rownames(total.All), cex=0.5, srt=90)
  dev.off()
  
}
