# 2014 analysis of all cases nationally

# start by extracting the necessary data (only human isolates)

library(dplyr)
library(lubridate)
library(xtable)

source("../../processing/helpers.R")

figure_dir <- "figures/epidemiology"
table_dir  <- "tables"

# read in data
data_file <- find_latest_version("../national_data")
data <- read.csv(data_file, stringsAsFactors=F)
data <- data %>% mutate(Date = as.Date(Date))
data <- data %>% mutate(YearMonth = (Year-2005)*12 + as.numeric(month(Date)))

pdf(file.path(figure_dir, "national_by_month.pdf"), width=8, height=4)
plot(Cases ~ YearMonth, type="l", data=data, col="black", xlim=range(YearMonth), xlab="", ylab="Cases per month", xaxs="i", xaxt="n")
#with(subset(cases_per_month, UR_bool == "Rural"), lines(YearMonth, value, col="red"))
# add the axis
years <- max(data$YearMonth)/12
axis(1, at=seq(0, years*12, by=12), label=rep("", years+1))
mtext(2005:2014, side=1, line=0.5, at=seq(6, years*12-6, 12))
#legend("topright", legend=c("Urban", "Rural"), col=1:2, lty="solid")
dev.off()

