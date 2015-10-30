# stuff for doing pictures for SA talks

# we first want a distribution of human STs and source STs

library(dplyr)
library(stringr)
library(tidyr)

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
