# stuff for doing pictures for SA talks

# we first want a distribution of human STs and source STs

library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)

dat = read.csv("~/data/R/islandR/running/20150603.csv", stringsAsFactors=F)

source_file <- "4_source_imputed"

sources <- read.csv(file.path("~/data/R/islandR/running", paste0(source_file, ".csv")), colClasses="character")
source_map <- as.numeric(sources$Number)
names(source_map) <- sources$DataSource

source_label_map <- unique(sources %>% dplyr::select(Number, Label))
source_labels <- str_replace(source_label_map$Label, "\\\\n", "\n")
names(source_labels) <- source_label_map$Number

# combine sources
dat = dat %>% mutate(SourceNumber = as.character(source_map[Source]))

dat = dat %>% left_join(source_label_map, by = c('SourceNumber' = 'Number'))
dat = dat %>% rename(ModelSource = Label) %>% mutate(ModelSource=ifelse(is.na(ModelSource), "Human", ModelSource))

# right, now generate some summaries of different types
test = dat %>% mutate(temp = 1L)

sts = dat %>% group_by(ST, ModelSource) %>% summarise(n=n()) %>% spread(ModelSource, n, fill=0)

write.csv(sts, "data/st_counts.csv", row.names=FALSE)

ur_sts = dat %>% filter(Source == "Human", !is.na(UR_bool)) %>% group_by(ST, UR_bool) %>% summarise(n=n()) %>% spread(UR_bool, n, fill=0)
write.csv(ur_sts, "data/ur_counts.csv", row.names=FALSE)

humans = dat %>% filter(Source == "Human") %>% mutate(Sampled.Date = as.Date(Sampled.Date), Month = month(Sampled.Date))

hbym = humans %>% group_by(Month, Year) %>% summarise(Count=n()) %>% mutate(Date = dmy(paste("01",Month,Year,sep="-"))) %>% ungroup %>% arrange(Date)

write.csv(hbym, "data/human_by_month.csv", row.names = FALSE)

# read in map stuff
library(splancs)
library(maptools)
library(RColorBrewer)
library(sp)
library(spatstat)
library(classInt)
library(dplyr)
darken <- function(col, d=0.5)
{
rgb(t(col2rgb(col)/255)*d, alpha=col2rgb(col, alpha=TRUE)[4]/255)
}
#' Input information
shape_files_folder <- "shape" # folder to load shapefiles from
continuous  <- FALSE      # set TRUE if the variable is continuous
num_colours <- 9         # shouldn't be too high
opacity     <- 0.8        # use to lighten the colours up a bit (1 = opaque, 0 = transparent)
shape_files_folder
shape_files_folder = "~/data/R/campy_db/shape"
map.shp <- readShapeSpatial(file.path(shape_files_folder, "midcentral_phu.shp"))
plot(map.shp)

u = unionSpatialPolygons(map.shp, rep(1, length(map.shp)))
writeSpatialShape(us, "test.shp")


# attribution mapping
attribution = read.csv("attribution.csv")
attribution = attribution %>% mutate(Year = (Time - 1) %/% 12 + 2005, Month = as.factor(Month))
attribution = attribution %>% left_join(humans %>% select(-Date))
write.csv(attribution, "data/attribution2.csv", row.names=FALSE)

attribution = read.csv("data/attribution.csv")
par(mfrow=c(4,1), mar=c(4,4,2,2))
for (src in levels(attribution$Source)) {
  d = attribution[attribution$Source == src,]
  n_times = max(d$Time)
  times = 1:n_times
  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(d$Source[1], "Urban", sep=" - "))
  polygon(c(times, rev(times)), c(d$ui[times], rev(d$li[times])), col="grey80", border=NA)
  lines(times, d$mu[times], lwd=2)

  plot(NULL, xlim=range(times), ylim=c(0,1), type="n", main=paste(d$Source[1], "Rural", sep=" - "))
  polygon(c(times, rev(times)), c(d$ui[times+n_times], rev(d$li[times+n_times])), col="grey80", border=NA)
  lines(times, d$mu[times+n_times], lwd=2)
}

par(mfrow=c(4,1), mar=c(4,4,2,2))
for (src in levels(attribution$Source)) {
  d = attribution[attribution$Source == src,]
  mu = d$mu * d$Count
  li = d$li * d$Count
  ui = d$ui * d$Count
  n_times = max(d$Time)
  times = 1:n_times
  plot(NULL, xlim=range(times), ylim=c(0,40), type="n", main=paste(d$Source[1], "Urban", sep=" - "))
  polygon(c(times, rev(times)), c(ui[times], rev(li[times])), col="grey80", border=NA)
  lines(times, mu[times], lwd=2)

  plot(NULL, xlim=range(times), ylim=c(0,40), type="n", main=paste(d$Source[1], "Rural", sep=" - "))
  polygon(c(times, rev(times)), c(ui[times+n_times], rev(li[times+n_times])), col="grey80", border=NA)
  lines(times, mu[times+n_times], lwd=2)
}

# TODO: Convert to dygraph?

# TODO: Map using leaflet?
