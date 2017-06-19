# spatial stuff
library(RColorBrewer)
#library(classInt)
library(maptools)

brewerBrBG9 <- rev(c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A"))

breaks <- c(-Inf, 0.2, 0.5, 0.8, 0.95, 1.05, 1.3, 2.0, 5.0, Inf)
col    = c("steelblue", "plum4")

get_col <- function(c, pal) {
  func <- colorRamp(pal, space="Lab")
  rgb(func(c), maxColorValue = 255)
}

ax_col <- "grey20"

temp_dir <- 'figures/epi'
dir.create(temp_dir)

load('data/epiclustR/TA2006.Rdata'); TA <- new

scases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases')
dim(scases) <- c(dim(scases)[1], prod(dim(scases)[2:3]))
weeks = as.Date(rownames(TA$data$cases))

# generate some interpolated random normals
set.seed(3)
d = sample(1:ncol(scases), 20)
d = c(d, d[1])
dim(scases[,d])

o = apply(scases[,d], 1, function(x) { spline(seq_along(x), x, method="periodic", xout=seq(1, length(x), length.out=length(x)*12+1))$y[-1]})

med <- apply(scases, 1, median)
lci <- apply(scases, 1, quantile, 0.25)
uci <- apply(scases, 1, quantile, 0.75)

for (i in seq_len(nrow(o))) {
  cat("running image", i, "of", nrow(o), "\n")
  png(file.path(temp_dir, sprintf("temp_fit%04d.png", i)), width=960, height=480, bg="transparent")
  par(mar=c(3,3,0.5,1), mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, apply(TA$data$cases, 1, sum), ylim=c(0,20), type='l', col="grey60", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='Cases', axes=FALSE, col.lab=ax_col)
  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, o[i,], lwd=2.5, col=alpha(col[1],1))
  lines(weeks, med, col=col[1], lwd=4)
  polygon(c(weeks, rev(weeks)), c(lci, rev(uci)), col=alpha(col[1], 0.5), border=NA)
  dev.off()
}


system(paste0("~/ffmpeg -y -r 24 -i ", file.path(temp_dir, 'temp_fit%04d.png'), " figures/temporal_fit.mp4"))

# spatial uncertainty of posterior for pre/post 2008
load('data/epiclustR/AU2006.Rdata'); AU <- new

U <- epiclustR:::ssapply(AU$mod, epiclustR:::extract_spatial, data=AU$data)
dim(U) <- c(dim(U)[1:2], prod(dim(U)[3:4]))

library(maptools)
phu <- readShapeSpatial('../nzsa_2016/maps/midcentral_phu')

num_cols <- 101
up_bks = quantile(U[U > 0], seq(1,num_cols,by=2)/num_cols)
lo_bks = quantile(U[U < 0], seq(num_cols-1,0,by=-2)/num_cols)
bks = round(pmax(up_bks, abs(lo_bks)), 2)
breaks = c(-rev(bks), bks)

alpha <- function(col, x = 0.5) {
  rgb(t(col2rgb(col)), alpha = x * 255, maxColorValue = 255)
}
brewerBrBG9 <- c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A")
#cols <- alpha(brewerBrBG9, 0.7)
cols <- colorRampPalette(brewerBrBG9, space='Lab')(num_cols)

num_alpha <- 10
breaks_alpha <- seq(0,1,length.out=num_alpha+1)

uncertain_func <- function(x) {
  y <- 1-(max(sum(x > 0), sum(x < 0)) / length(x)*2 - 1)
  # those more uncertain, make more so (y == 1)
  sqrt(y)
}

uncertainty <- apply(U, 1:2, uncertain_func)
uncertainty <- apply(uncertainty, 2, function(x) { as.numeric(cut(x, breaks=breaks_alpha, include.lowest = TRUE)) })

cols_unc <- lapply(cols, function(x) { colorRampPalette(c(x, 'grey50'), space='Lab')(num_alpha+1) })
alpha_cols <- simplify2array(cols_unc)

if (is.null(bbox) || !is.matrix(bbox))
  bbox = sp::bbox(phu)

plot_map <- function(u, uncertainty, bbox = sp::bbox(phu)) {
  spat_risk <- cbind(AU$data$spat_list, Risk=u)
  spat_risk$Spatial = as.numeric(as.character(spat_risk$Spatial))
  
  map_dat <- phu@data %>%
    dplyr::left_join(spat_risk, by=c('MB06' = 'Spatial'))
  
  col_val <- cut(map_dat$Risk, breaks = breaks)
#  print(head(cbind(uncertainty, col_val)))
#  print(dim(cbind(uncertainty, col_val)))
#  print(dim(alpha_cols))

  map_col <- alpha(alpha_cols[cbind(uncertainty, col_val)], 1)
  
  sp::plot(phu, col=map_col, lwd=0.02, border='grey80', xlim=bbox[1,], ylim=bbox[2,])
}

set.seed(5)
iters <- sample(seq_along(U[1,1,]), 40)
V <- apply(U[,,iters], 1:2, function(x) { spline(seq_len(length(x)+1), c(x, x[1]), method='periodic', xout=seq(1, length(x)+1, by=1/20)[-1])$y })
for (i in seq_len(dim(V)[1])) {
  png(file.path(temp_dir, sprintf("spatial_fit%04d.png", i)), width=960, height=480)
  par(mfrow=c(1,2), mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", dim(V)[1], "\n")
  lapply(seq_len(dim(V)[3]), function(y) { plot_map(V[i,,y], uncertainty[,y]) })
  dev.off()
}
#system("convert -delay 4 -loop 0 -dispose background spatial_fit*.png spatial_fit2.gif")
#system("convert spatial_fit2.gif -transparent white spatial_fit.gif")

system(paste0("~/ffmpeg -y -r 24 -i ", file.path(temp_dir, 'spatial_fit%04d.png'), " figures/spatial_fit.mp4"))
system(paste0("avconv -y -r 24 -i ", file.path(temp_dir, 'spatial_fit%04d.png'), " figures/spatial_fit.mp4"))

#system("avconv -y -r 24 -i spatial_fit%04d.png spatial_fit.mp4")

for (i in seq_len(dim(V)[1])) {
  png(file.path(temp_dir, sprintf("spatial_palmy_fit%04d.png", i)), width=960, height=480)
  par(mfrow=c(1,2), mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", dim(V)[1], "\n")
  lapply(seq_len(dim(V)[3]), function(y) { plot_map(V[i,,y], uncertainty[,y], bbox=matrix(c(2727946,6086900,2734889,6094322), 2)) })
  dev.off()
}
#system("convert -delay 4 -loop 0 -dispose background spatial_palmy_fit*.png spatial_palmy_fit2.gif")
#system("convert spatial_palmy_fit2.gif -transparent white spatial_palmy_fit.gif")

system(paste0("~/ffmpeg -y -r 24 -i ", file.path(temp_dir, 'spatial_palmy_fit%04d.png'), " figures/spatial_fit_palmy.mp4"))
system(paste0("avconv -y -r 24 -i ", file.path(temp_dir, 'spatial_palmy_fit%04d.png'), " figures/spatial_fit_palmy.mp4"))

#system("avconv -y -r 24 -i spatial_palmy_fit%04d.png spatial_palmy_fit.mp4")

plot_temporal <- function(weeks, scases, ecases) {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, ecases, ylim=c(0,12), type='l', col="red", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, scases, col="black", lwd=2)
#  rect(weeks[week]-14, 0, weeks[week]+14, 12, col=alpha("steelblue", 0.7), border=NA)
}

scasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
scasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
weeks = as.Date(rownames(TA$data$cases))

plot_temporal(weeks, scasesTA, ecasesTA)
plot_temporal(weeks, scasesAU, ecasesAU)

# TODO: Merge the AU in with the TA

setdiff(which(ecasesAU - scasesAU > 1), which(ecasesTA - scasesTA > 1))




# grab outreak data
roll_outbreak_probs <- function(mod, window=20) {
  X <- epiclustR:::ssapply(mod, epiclustR:::extract_variable, 'X')
  mX = apply(X, 1:2, mean)
  
  # make it a rolling window
  mR <- mX
  for (i in 1:window)
    mR[-(1:i),] <- mR[-(1:i),] + mX[1:(nrow(mX)-i),] * exp(-5/window*i)
  mR
}
Rta <- roll_outbreak_probs(TA$mod)
Rau <- roll_outbreak_probs(AU$mod)


# plot these through time
plot_ob_map <- function(sl, u, x, bbox = sp::bbox(phu)) {
  spat_risk <- cbind(sl, Risk=u)
  outbreaks <- data.frame(Region = 1:length(x), P=x)

  map_dat <- phu@data %>%
    dplyr::left_join(spat_risk, by=c('MB06' = 'Spatial')) %>%
    dplyr::left_join(outbreaks, by="Region")

  alpha2 <- function(col, x = 0.5) {
    rgb(t(col2rgb(col))*x + (1-x) * 255, maxColorValue = 255)
  }

  vals <- cut(map_dat$Risk, breaks = breaks)
  map_col <- alpha2(cols[vals], 0.7)

  red <- function(x) {
    red_func <- colorRamp(brewer.pal(9, "Reds"), space="Lab")
    rgb(red_func(pmin(x,1)), maxColorValue = 255)
  }
  sp::plot(phu, col=red(map_dat$P), lwd=0.02, border='grey80', xlim=bbox[1,], ylim=bbox[2,])
#  sp::plot(phu, add=TRUE, col=red(map_dat$P), border=NA)
}


# plot these through time
plot_ob_map2 <- function(map, x, bbox = sp::bbox(phu)) {
  red <- function(x) {
    red_func <- colorRamp(brewer.pal(9, "Reds")[-1], space="Lab")
    rgb(red_func(pmin(x,1)), maxColorValue = 255)
  }
  sp::plot(map, col=red(x), lwd=0.5, border='grey50', xlim=bbox[1,], ylim=bbox[2,])
  #  sp::plot(phu, add=TRUE, col=red(map_dat$P), border=NA)
}

plot_temporal <- function(weeks, scases, ecases, week) {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, ecases, ylim=c(0,12), type='l', col="red", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, scases, col="black", lwd=2)
  rect(weeks[week]-14, 0, weeks[week]+14, 12, col=alpha("steelblue", 0.7), border=NA)
}

library(maptools)
phu <- readShapeSpatial('maps/midcentral_phu')
AU$data$spat_list$Region <- as.numeric(as.factor(AU$data$spat_list$Region))

num_cols <- 101
up_bks = quantile(U[U > 0], seq(1,num_cols,by=2)/num_cols)
lo_bks = quantile(U[U < 0], seq(num_cols-1,0,by=-2)/num_cols)
bks = round(pmax(up_bks, abs(lo_bks)), 3)
breaks = c(-rev(bks), bks)

au_reg = phu@data %>% left_join(AU$data$spat_list, by=c("MB06"="Spatial"))
au_shp = unionSpatialPolygons(phu, sprintf("%02d", au_reg$Region))
ta_reg = phu@data %>% left_join(TA$data$spat_list, by=c("MB06"="Spatial"))
ta_shp = unionSpatialPolygons(phu, sprintf("%02d", ta_reg$Region))

scasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
scasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'ecases'), 1, median)

library(dplyr)

for (i in 1:nrow(Rta)) {#seq_len(nrow(U))) {
  png(sprintf("png/outbreak_map_nt_%04d.png", i), width=960, height=480, bg="transparent")
  layout(matrix(1:4, 2, 2), heights=c(4,1), widths=c(1,1))
  par(mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", nrow(Rta), "\n")
  # tail mX off over a bunch of frames
  plot_ob_map2(ta_shp, Rta[i,])
  par(mai=c(0.5,0.416,0,0.416))
  plot_temporal(weeks, scasesTA, ecasesTA, i)
  par(mai=c(0.1,0.5,0.1,0.5))
  plot_ob_map2(au_shp, Rau[i,], bbox=matrix(c(2717946,6078900,2754889,6102322), 2))
  par(mai=c(0.5,0.416,0,0.416))
  plot_temporal(weeks, scasesAU, ecasesAU, i)
  # plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
  dev.off()
}

system("~/ffmpeg -y -r 10 -i png/outbreak_map_nt_%04d.png -b:v 1M video/outbreak_map.mp4")

system("convert -delay 10 -loop 0 -dispose background outbreak_map_nt_*.png outbreak_map2.gif")
system("convert outbreak_map2.gif -transparent white outbreak_map.gif")

animation::saveVideo(
for (i in 1:nrow(Uta)) {#seq_len(nrow(U))) {
#  png(sprintf("outbreak_map_nt_%04d.png", i), width=960, height=480, bg="transparent")
  layout(matrix(1:4, 2, 2), heights=c(4,1), widths=c(1,1))
  par(mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", nrow(U), "\n")
  # tail mX off over a bunch of frames
  plot_ob_map(TA$data$spat_list, Uta[i,], Rta[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesTA, ecasesTA, i)
  par(mai=c(0,0,0,0))
  plot_ob_map(AU$data$spat_list, Uau[i,], Rau[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesAU, ecasesAU, i)
  # plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
#  dev.off()
}, video.name="outbreak_map.mp4", ffmpeg="avconv", ani.width=960, ani.height=480, interval=1/24)


system("convert -delay 10 -loop 0 -dispose background outbreak_map_nt_*.png outbreak_map2.gif")
system("convert outbreak_map2.gif -transparent white outbreak_map.gif")


# TODO: could do just one, combined plot using orange + purple instead of red
i <- 137
png(sprintf("test.png", i), width=960, height=480, bg="transparent")
layout(matrix(c(1,3,2,3), 2, 2), heights=c(4,1))
par(mai=c(0,0,0,0), bg="#FFFFFF")
cat("plotting", i, "of", nrow(U), "\n")
# tail mX off over a bunch of frames
plot_ob_map(TA$data$spat_list, Uta[i,], Rta[i,])
plot_ob_map(AU$data$spat_list, Uau[i,], Rau[i,])
par(mai=c(0.5,0.5,0,0.5))
plot_temporal2(weeks, scasesTA, ecasesTA, ecasesAU, i)
# plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
dev.off()

plot_temporal2 <- function(weeks, scases, ecases1, ecases2, week) {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  PuOr <- brewer.pal(11, "PuOr")[c(2,10)]
  plot(weeks, ecases1, ylim=c(0,12), type='n', col=PuOr[1], xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  segments(weeks, y0=scases, y1=ecases1, col=PuOr[1], lwd=2)
  segments(weeks, y0=scases, y1=ecases2, col=PuOr[2], lwd=2)
#  lines(weeks, ecases2, col=puor[2], lwd=1)
  lines(weeks, scases, col="black", lwd=2)
  rect(weeks[week]-15, 0, weeks[week]+15, 12, col=alpha("steelblue", 0.7), border=NA)
}

# TODO: Add the temporal plot underneath this...
scasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
scasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'ecases'), 1, median)

plot_temporal(weeks, scasesTA, ecasesTA)
plot_temporal(weeks, scasesAU, ecasesAU)

plot_temporal(weeks, scasesTA, ecasesTA)
axis(2)

# Find the Raw Milk outbreak and Pahiatua outbreaks
wch_rawmilk <- which(ecasesTA - scasesTA > 0.4)[5]
wch_pahiatua <- which(ecasesTA - scasesTA > 0.4)[1]

ecasesTA[wch_rawmilk] <- scasesTA[wch_rawmilk]+1

# make a plot with these on them
cases <- apply(TA$data$cases, 1, sum)

png("figures/14_outbreak_plot.png", width=fig_width, height=fig_height)
par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
plot(weeks, cases, ylim=c(0,20), type='l', col="grey60", xaxs='i', yaxs='i', xlab='', lwd=1,
     ylab='', axes=FALSE, col.lab=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
lines(weeks, ecasesTA, col="red", lwd=1.5)
lines(weeks, scasesTA, col="black", lwd=3)
dev.off()

png("figures/14_outbreak_plot_highlight.png", width=fig_width, height=fig_height)
par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
plot(weeks, cases, ylim=c(0,20), type='l', col="grey60", xaxs='i', yaxs='i', xlab='', lwd=1,
     ylab='', axes=FALSE, col.lab=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
lines(weeks, ecasesTA, col="red", lwd=1.5)
lines(weeks, scasesTA, col="black", lwd=3)
dim <- 0.85
bright <- 0.3
rect(weeks[1],0,weeks[wch_pahiatua-5],20,col=alpha("white", dim), border=NA)
rect(weeks[wch_pahiatua-5],0,weeks[wch_pahiatua+5],20,col=NA, border='black', lwd=3)
rect(weeks[wch_pahiatua+5],0,weeks[wch_rawmilk-5],20,col=alpha("white", dim), border=NA)
rect(weeks[wch_rawmilk-5],0,weeks[wch_rawmilk+5],20,col=NA, border='black', lwd=3)
rect(weeks[wch_rawmilk+5],0,weeks[length(weeks)],20,col=alpha("white", dim), border=NA)
dev.off()



rect(weeks[week]-14, 0, weeks[week]+14, 12, col=alpha("steelblue", 0.7), border=NA)


TA$data$cases[wch_rawmilk,]

TA$data$mbrg
case_list_ta <- TA$data$case_list %>% mutate(Spatial = as.character(Spatial)) %>% left_join(TA$data$spat_list)


rawmilk_cases <- case_list_ta %>% filter(ReportWeek %in% weeks[wch_rawmilk+-1:1], Region == "Palmerston North City")

# join up to main database to get info on ST
library(campydb)
s <- read_samples("~/data/db_config.txt") %>% left_join(read_isolates("~/data/db_config.txt"), by="LabID")
e1 <- read_episurv_extra("~/data/db_config.txt")
db <- s %>% left_join(e1, by=c("HospitalNo" = "hospitalnu"))

rawmilk_final <- rawmilk_cases %>% left_join(db, by=c("CaseID" = "episurvnum")) %>%
  select(one_of(names(rawmilk_cases)), ST) %>% filter(!is.na(ST))

write.csv(rawmilk_final, "data/outbreak_rawmilk.csv", row.names=FALSE)



# compute the expected cases for all observations
scases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases')
ecases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases')
dim(scases) <- c(dim(scases)[1], prod(dim(scases)[2:3]))
dim(ecases) <- c(dim(ecases)[1], prod(dim(ecases)[2:3]))
weeks = as.Date(rownames(TA$data$cases))

# generate some interpolated random normals
set.seed(3)
d = sample(1:ncol(scases), 20)
d = c(d, d[1])

o_s = apply(scases[,d], 1, function(x) { spline(seq_along(x), x, method="periodic", xout=seq(1, length(x), length.out=length(x)*12+1))$y[-1]})
o_e = apply(ecases[,d], 1, function(x) { spline(seq_along(x), x, method="periodic", xout=seq(1, length(x), length.out=length(x)*12+1))$y[-1]})

med <- apply(scases, 1, median)
lci <- apply(scases, 1, quantile, 0.25)
uci <- apply(scases, 1, quantile, 0.75)

med_e <- apply(ecases, 1, median)

for (i in seq_len(nrow(o))) {
  png(sprintf("ob_fit%04d.png", i), width=960, height=480)
  par(mar=c(3,3,0.5,1), mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, apply(TA$data$cases, 1, sum), ylim=c(0,20), type='l', col="grey70", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='Cases', axes=FALSE, col.lab=ax_col)
  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, o_e[i,], lwd=1, col=alpha("red",0.5))
  lines(weeks, o_s[i,], lwd=1, col=alpha("black",1))
  lines(weeks, med_e, col="red", lwd=1)
  lines(weeks, med, col="black", lwd=2)
#  polygon(c(weeks, rev(weeks)), c(lci, rev(uci)), col=alpha(col[1], 0.3), border=NA)
  dev.off()
}

system("convert -delay 4 -loop 0 -dispose background ob_fit*.png outbreak_fit2.gif")
system("convert outbreak_fit2.gif -transparent white outbreak_fit.gif")



###########################

cases <- apply(TA$data$cases, 1, sum)
dates <- as.Date(names(cases))
mod <- loess(cases ~ as.numeric(dates), span = 1/10)
y <- predict(mod, se=TRUE)

png("figures/11_temporal_cases.png", width = fig_width, height = fig_height)
par(mar=c(3,3,0.5,1), mgp=c(2,.7,0), tck=-.015)
plot(cases ~ dates, type="l", col="grey60", axes=FALSE, xaxs="i", yaxs="i", ylim=c(0,20), xlab="", ylab="Cases", col.lab=ax_col)
axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
y_fit = y$fit
polygon(c(dates, rev(dates)), c(y_fit+y$se.fit, rev(y_fit - y$se.fit)), col=alpha(col[1], 0.5), border=NA)
lines(dates, y_fit, col=col[1], lwd=2)
dev.off()


############ SPATIAL CASES

cases1 <- apply(TA$data$cases[TA$data$t2p==1,], 2, sum)/sum(TA$data$t2p==1)*52
cases2 <- apply(TA$data$cases[TA$data$t2p==2,], 2, sum)/sum(TA$data$t2p==2)*52

spat_case <- cbind(TA$data$spat_list, Case1 = cases1 / colSums(TA$data$popn), Case2 = cases2 / colSums(TA$data$popn))
map_data <- phu@data %>% mutate(MB06 = as.character(MB06)) %>% dplyr::left_join(spat_case, by=c('MB06' = 'Spatial')) %>%
  mutate(Case1 = ifelse(is.na(Case1), 0, Case1), Case2 = ifelse(is.na(Case2), 0, Case2))

# do quantile stuff
brewerBrBG9 <- c("#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A")
rate <- c(map_data$Case1,map_data$Case2)
cuts <- quantile(log(rate[rate > 0]), 0:5/5)
cols1 <- brewerBrBG9[cut(log(map_data$Case1), cuts)]
cols1[is.na(cols1)] <- brewerBrBG9[1]
cols2 <- brewerBrBG9[cut(log(map_data$Case2), cuts)]
cols2[is.na(cols2)] <- brewerBrBG9[1]

png("figures/12_spatial_cases08.png", width = fig_width, height = fig_height)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(phu, col=cols1, border="grey60", lwd=0.3)
plot(phu, xlim=c(2727946, 2734889), ylim=c(6086900, 6094322), col=cols1, border="grey60", lwd=0.3)
dev.off()

png("figures/13_spatial_cases16.png", width = fig_width, height = fig_height)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(phu, col=cols2, border="grey60", lwd=0.3)
plot(phu, xlim=c(2727946, 2734889), ylim=c(6086900, 6094322), col=cols2, border="grey60", lwd=0.3)
dev.off()


