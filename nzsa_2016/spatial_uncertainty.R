# spatial stuff
library(RColorBrewer)
library(classInt)
library(maptools)

brewerBrBG9 <- rev(c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A"))

breaks <- c(-Inf, 0.2, 0.5, 0.8, 0.95, 1.05, 1.3, 2.0, 5.0, Inf)

shapeFile <- "~/Massey/Projects/EpiclustR/PrepareDatasets/NZ/Data/NZ_shapefiles/AU2006.shp"
map.shp <- readShapeSpatial(shapeFile)
ids <- slot(map.shp, "data")[[1]]
map.shp$num <- 1:length(ids);

xylims <- attr(map.shp, "bbox")
xylims2 <- matrix(c(1752500, 1757500, 5918500, 5922000), 2, 2, byrow=T)
xylims <- xylims2 + c(-3000,-3000,3000,3000)

num <- c(229:235, 245:257, 284:285, 293)

xylims <- bbox(map.shp[num,])

ratio <- (xylims[2,2] - xylims[2,1]) / (xylims[1,2] - xylims[1,1])

get_col<-function(c, pal)
{
  o <- ((floor(c) + 1 - c)*col2rgb(pal[floor(c)]) + (c-floor(c))*col2rgb(pal[floor(c)+1]))/255
  return(rgb(o[1,], o[2,], o[3,]))
}

set.seed(2)
risk <- runif(6,1,4)

# generate a bunch of random risks around the average risk
set.seed(3)
d = rnorm(20, mean=mean(risk), sd=1)
d = c(d, d[1])
dt = c(0, cumsum(abs(diff(d))))
o = spline(dt, d, method="periodic", xout = seq(min(dt), max(dt), length.out=length(d)*12+1))$y[-1]

for (i in seq_along(o)) {
  png(sprintf("plot%04d.png", i), width=640, height=400, bg="transparent")
  par(pin = c(5/sqrt(ratio), sqrt(ratio) * 5), mar=c(0,0,0,0), bg="#FFFFFF")
  plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0, bty="n")
  
  cols <- rep("grey90",length(ids))
  cols[251] <- get_col(o[i], brewerBrBG9)
  cols[c(248:250,252:253,229)] <- get_col(risk, brewerBrBG9)

  plot(map.shp[num,], lty=1, col=cols[num], border="black", lwd=0.6, add=T)
  
  coords <- slot(slot(map.shp, "polygons")[[251]], "labpt")
  text(coords[1], coords[2], expression(U[i]))
  dev.off()
}

system("convert -delay 4 -loop 0 plot*.png spatial2.gif")
system("convert spatial2.gif -transparent white spatial.gif")



# temporal uncertainty here...
t = c(-2,-1,1,2)
R = c(3,1.5,2,5)
m = lm(R ~ poly(t, 2))
that = seq(-2,2,by=0.02)
Rt = predict(m, data.frame(t=0))
quad_fit = predict(m, data.frame(t=that))
y = seq(-3,3,by=0.03)
# generate some interpolated random normals
set.seed(3)
d = rnorm(20, sd=0.3)
d = c(d, d[1])
dt = c(0, cumsum(abs(diff(d))))
o = spline(dt, d, method="periodic", xout = seq(min(dt), max(dt), length.out=length(d)*12+1))$y[-1]
foo = animation::saveGIF(for (i in o) {
  plot(R ~ t, axes=FALSE, bty='n', ylab='', xlab='', ylim=c(0,5), cex=1.3)
  polygon(c(dnorm(y)*0.3, -dnorm(y)*0.3), Rt + 0.3*c(y, rev(y)), col="grey70", border=NA)
  lines(quad_fit ~ that, lwd=2)
  points(0, Rt, pch=19, cex=1.3)
  axis(1, at=-2:2, labels = c(expression(t-2, t-1, t, t+1, t+2)))
  tr <- -2:2
  Rr <- predict(m, data.frame(t = tr))
  Rr[3] = Rr[3] + i
  mr <- lm(Rr ~ poly(tr, 4))
  lines(predict(mr, data.frame(tr=that)) ~ that, col="red", lwd=2)
  points(0, Rr[3], pch=19, cex=1.3, col="red")
}, interval=1/20, movie.name="temporal.gif", ani.width=640, ani.height=320)



ax_col <- "grey20"

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
  png(sprintf("temp_fit%04d.png", i), width=960, height=480, bg="transparent")
  par(mar=c(3,3,0.5,1), mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, apply(TA$data$cases, 1, sum), ylim=c(0,20), type='l', col="grey60", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='Cases', axes=FALSE, col.lab=ax_col)
  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, o[i,], lwd=1.5, col=alpha(col[1],1))
  lines(weeks, med, col=col[1], lwd=2.5)
  polygon(c(weeks, rev(weeks)), c(lci, rev(uci)), col=alpha(col[1], 0.5), border=NA)
  dev.off()
}

system("convert -delay 4 -loop 0 -dispose background temp_fit*.png temporal_fit2.gif")
system("convert temporal_fit2.gif -transparent white temporal_fit.gif")




# spatial uncertainty of posterior for pre/post 2008
U <- epiclustR:::ssapply(TA$mod, epiclustR:::extract_spatial, data=TA$data)
dim(U) <- c(dim(U)[1:2], prod(dim(U)[3:4]))

library(maptools)
phu <- readShapeSpatial('maps/midcentral_phu')

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
cols <- colorRampPalette(brewerBrBG9)(num_cols)

if (is.null(bbox) || !is.matrix(bbox))
  bbox = sp::bbox(phu)

plot_map <- function(u, bbox = sp::bbox(phu)) {
  spat_risk <- cbind(TA$data$spat_list, Risk=u)
  
  map_dat <- phu@data %>%
    dplyr::left_join(spat_risk, by=c('MB06' = 'Spatial'))
  
  vals <- cut(map_dat$Risk, breaks = breaks)
  map_col <- alpha(cols[vals], 1)
  
  sp::plot(phu, col=map_col, lwd=0.02, border='grey80', xlim=bbox[1,], ylim=bbox[2,])
}

set.seed(5)
iters <- sample(seq_along(U[1,1,]), 40)
V <- apply(U[,,iters], 1:2, function(x) { spline(seq_len(length(x)+1), c(x, x[1]), method='periodic', xout=seq(1, length(x)+1, by=1/20)[-1])$y })
for (i in seq_len(dim(V)[1])) {
  png(sprintf("spatial_fit%04d.png", i), width=960, height=480)
  par(mfrow=c(1,2), mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", dim(V)[1], "\n")
  apply(V[i,,], 2, function(y) { plot_map(y) })
  dev.off()
}
#system("convert -delay 4 -loop 0 -dispose background spatial_fit*.png spatial_fit2.gif")
#system("convert spatial_fit2.gif -transparent white spatial_fit.gif")
system("avconv -y -r 24 -i spatial_fit%04d.png spatial_fit.mp4")

for (i in seq_len(dim(V)[1])) {
  png(sprintf("spatial_palmy_fit%04d.png", i), width=960, height=480)
  par(mfrow=c(1,2), mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", dim(V)[1], "\n")
  apply(V[i,,], 2, function(y) { plot_map(y, bbox=matrix(c(2727946,6086900,2734889,6094322), 2)) })
  dev.off()
}
#system("convert -delay 4 -loop 0 -dispose background spatial_palmy_fit*.png spatial_palmy_fit2.gif")
#system("convert spatial_palmy_fit2.gif -transparent white spatial_palmy_fit.gif")

system("avconv -y -r 24 -i spatial_palmy_fit%04d.png spatial_palmy_fit.mp4")




# outbreak plot...
mean_case_rate <- function(mod, data) {
  ecases <- matrix(0, nrow(data$cases), ncol(data$cases))
  rownames(ecases) <- rownames(data$cases)
  colnames(ecases) <- colnames(data$cases)
  for (i in seq_along(mod)) {
    cat("up to chain", i, "\n")
    for (j in seq_along(mod[[i]])) {
      ecases <- ecases + epiclustR::log_case_rate(data, mod[[i]][[j]], smoothed=TRUE)
    }
  }
  ecases <- ecases / (length(mod) * length(mod[[1]]))
  ecases - mean(ecases)
}

Uta <- mean_case_rate(TA$mod, TA$data)
Uau <- mean_case_rate(AU$mod, AU$data)

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
  rect(weeks[week]-15, 0, weeks[week]+15, 12, col=alpha("steelblue", 0.7), border=NA)
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

library(dplyr)

for (i in 1:nrow(Uta)) {#seq_len(nrow(U))) {
  png(sprintf("outbreak_map_nt_%04d.png", i), width=960, height=480, bg="transparent")
  layout(matrix(1:4, 2, 2), heights=c(4,1), widths=c(1,1))
  par(mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", nrow(U), "\n")
  # tail mX off over a bunch of frames
  plot_ob_map2(ta_shp, Rta[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesTA, ecasesTA, i)
  par(mai=c(0,0,0,0))
  plot_ob_map2(au_shp, Rau[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesAU, ecasesAU, i)
  # plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
  dev.off()
}

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

system("avconv -y -r 10 -i outbreak_map_nt_%04d.png -b:v 1M outbreak_map.mp4")

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











pdf("priors_u1.pdf", width=5/sqrt(ratio), height=5*sqrt(ratio))
par(pin = c(4.6/sqrt(ratio), sqrt(ratio) * 4.6), omi = c(0,0,0,0))
plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0)

cols <- rep("grey90",length(ids))
cols[251] <- "grey60"
cols[c(248:250,252:253,229)] <- get_col(risk, pal)
num <- c(229:235, 245:257, 284:285, 293)
plot(map.shp[num,], lty=1, col=cols[num], border="black", lwd=0.6, add=T)

coords <- slot(slot(map.shp, "polygons")[[251]], "labpt")
text(coords[1], coords[2], expression(U[i]))
box()
dev.off()

risk <- runif(6,1,4)

pdf("priors_u2.pdf", width=5/sqrt(ratio), height=5*sqrt(ratio))
par(pin = c(4.6/sqrt(ratio), sqrt(ratio) * 4.6), omi = c(0,0,0,0))
plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0)

cols <- rep("grey90",length(ids))
cols[251] <- get_col(mean(risk), pal)
cols[c(248:250,252:253,229)] <- get_col(risk, pal)
num <- c(229:235, 245:257, 284:285, 293)
plot(map.shp[num,], lty=1, col=cols[num], border="black", lwd=0.6, add=T)

coords <- slot(slot(map.shp, "polygons")[[251]], "labpt")
text(coords[1], coords[2], expression(U[i]))

box()
dev.off()

pdf("priors_u3.pdf", width=5/sqrt(ratio), height=5*sqrt(ratio))
par(pin = c(4.6/sqrt(ratio), sqrt(ratio) * 4.6), omi = c(0,0,0,0))
plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0)

cols <- rep("grey90",length(ids))
cols[251] <- get_col(mean(risk)-1, pal)
cols[c(248:250,252:253,229)] <- get_col(risk, pal)
num <- c(229:235, 245:257, 284:285, 293)
plot(map.shp[num,], lty=1, col=cols[num], border="black", lwd=0.6, add=T)

coords <- slot(slot(map.shp, "polygons")[[251]], "labpt")
text(coords[1], coords[2], expression(U[i]))
box()
dev.off()

pdf("priors_u4.pdf", width=5/sqrt(ratio), height=5*sqrt(ratio))
par(pin = c(4.6/sqrt(ratio), sqrt(ratio) * 4.6), omi = c(0,0,0,0))
plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0)

cols <- rep("grey90",length(ids))
cols[251] <- get_col(mean(risk)+1, pal)
cols[c(248:250,252:253,229)] <- get_col(risk, pal)
num <- c(229:235, 245:257, 284:285, 293)
plot(map.shp[num,], lty=1, col=cols[num], border="black", lwd=0.6, add=T)

coords <- slot(slot(map.shp, "polygons")[[251]], "labpt")
text(coords[1], coords[2], expression(U[i]))
box()
dev.off()

zhome <- "~/Massey/Projects/EpiclustR/runs/NZ_fullrun/"
workingDir         <- paste(zhome,"NZ/RUX2_region",sep="") # directory containing this code
casesFile          <- paste(zhome,"NZ/Data.txt",sep="")
regionFile         <- paste(zhome,"NZ/Regions.txt",sep="")      # Region file (single row for each meshblock with region number)
regionNamesFile    <- paste(zhome,"NZ/RegionNames.txt",sep="")  # Region number:name:bigname mapping
meshblocksFile     <- paste(zhome,"NZ/Meshblocks.txt", sep="")
temporalFile       <- "smoothedCases.txt"
spatiotemporalFile <- "expectedCases.txt"
outputFile         <- "pnorth.pdf"                              # output file
startDate          <- "2008-08-30"                               # date of first timepoint

library(RColorBrewer)
setwd(workingDir)


# ignore the burnin samples
burn_in_samples <- 5000/20

# to get the latter we'd need a running cummulative "RU" similar to cumX.txt

expectedCasesRU <- read.table(temporalFile)
expectedCasesRU <- apply(expectedCasesRU[-(1:burn_in_samples),],2,mean)

num_times <- length(expectedCasesRU)

cases <- matrix(scan(casesFile), nrow=num_times)
observedCases <- apply(cases,1,sum)

expectedCasesRUX <- read.table(spatioTemporalFile)
expectedCasesRUX <- apply(expectedCasesRUX[-(1:burn_in_samples),],2,mean)

X <- as.matrix(read.table("posteriorX.txt"))
tolerance <- quantile(X,c(0.75,0.95,0.98,0.99))

regions <- scan(regionFile)
if (file.exists(regionNamesFile))
{
  region_names <- read.table(regionNamesFile, header=T)
} else
{
  region_names <- data.frame(number=1:num_regions, region=1:num_regions, bigregion=1:num_regions);
}
num_regions <- max(regions)

# we want the actual expected cases for each region at each time point for the last 4 weeks
# ideally this would be done by knowing fe+R_t+U_i for each region at each iteration (or the cummulative sum thereof)

# instead, we do it inefficiently...

mbs <- read.table(meshblocksFile)
num_mbs <- length(regions)

R <- matrix(scan("R.txt", skip=burn_in_samples), ncol=num_times, byrow=T)
fe <- scan("fixedEffects.txt", skip=burn_in_samples)
U <- matrix(scan("U.txt", skip=burn_in_samples), ncol=num_mbs, byrow=T)

wch <- list()
for (r in 1:num_regions)
{
  wch[[r]] <- which(regions==r)
}

num_weeks <- min(4,num_times)

expectedCases <- matrix(0, num_regions, num_weeks)
observedCases <- matrix(0, num_regions, num_weeks)
for (t in 1:num_weeks)
{
  e <- rep(0,length(num_mbs))
  for (m in 1:num_mbs)
  {
    e[m] <- mbs[m,2]*mean(exp(fe + R[,num_times-num_weeks+t] + U[,m]))
  }
  for (r in 1:num_regions)
  {
    expectedCases[r,t] <- sum(e[wch[[r]]])
    observedCases[r,t] <- sum(cases[num_times-num_weeks+t,wch[[r]]])
  }
}


# find where to place our year labels
year_labels <- NULL
for (t in 1:num_times)
{
  d <- as.POSIXlt(as.Date(startDate) + 7*(t-1))
  week <- floor(d$yday/7)
  if (week == 25)
    year_labels <- rbind(year_labels,c(t,d$year+1900))
}

# find the first year
first_week <- floor(as.POSIXlt(as.Date(startDate))$yday/7)
num_years <- (num_times+52) %/% 52


# now output each region


reg <- 1:num_regions

reg <- 1

reg <- 39

pdf("PN_outbreak2.pdf", width=8, height=5)
max_cases_per_region <- round(quantile(cases_per_region, 0.99) / 3 + 0.5)*3
for (r in reg)
{
  m <- max(max_cases_per_region, cases_per_region[,r])
  title <- paste(region_names$region[r],", ",region_names$bigregion[r], sep="")
  plot(1:num_times, t="n", xlab="", xaxt="n", yaxt="n", ylab="", main=title, xlim=c(1,num_times), ylim=c(-1,1))
  polygon(141.2+c(0,2.6,2.6,0),c(-2,-2,2,2),col="grey80", border="grey80")
  axis(2, labels=c(0,0.5,1), at=c(0,0.5,1))
  segments(1:num_times, 0, 1:num_times, as.numeric(X[r,]), col="black")
  case_axis <- seq(0,max_cases_per_region,length.out=4)
  axis(4, labels=case_axis, at=-1*case_axis/m)
  wh <- which(cases_per_region[,r] > 0)
  segments(wh, 0, wh, -cases_per_region[wh,r]/m, col="green3")

  axis(1, at=seq(-first_week,num_years*52-first_week,by=52), labels = FALSE)
  axis(1, at=seq(-first_week,num_years*52-first_week,by=52/12), labels = FALSE, tck=-.02)
  for (i in 1:nrow(year_labels))
    mtext(year_labels[i,2], side=1, line=1, at=year_labels[i,1], cex=0.9)
  legend(x = "topleft", fill= c("black", "green3"), legend=c("Outbreak probability", "Observed cases"), bg = "white", bty="n")
}
dev.off()


