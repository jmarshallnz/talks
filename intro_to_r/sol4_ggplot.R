#' ---
#' title: Intro to R - Lab 4 using ggplot2
#' author: Jonathan
#' ---
#' 
#' ## Exercise 1
#' 
#' Read in the data
goats = read.csv("http://www.massey.ac.nz/~jcmarsha/rcourse/goatsdata1.csv")

#' Load the ggplot2 library. If you haven't already, install it via the packages menu
#' in RStudio (bottom right)
library(ggplot2)

#' Now plot - just the data first
ggplot(goats, aes(x=chest, y=weight)) +
  geom_point()

#' Use a different plotting symbol (shape) and change the axis labels
ggplot(goats, aes(x=chest, y=weight)) +
  geom_point(shape=3) +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' Colour by gender. Note that because gender is numeric, it will
#' produce a colour scale rather than a discrete colour palette
ggplot(goats, aes(x=chest, y=weight, col=gender)) +
  geom_point() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' We can switch gender to a `factor` to get a discrete colour
#' palette
ggplot(goats, aes(x=chest, y=weight, col=factor(gender))) +
  geom_point() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' Modifying the legend is done by using `scale_color_discrete`
#' This is because it is a scale (an aesthetic always has a scale)
#' and that scale is for color (rather than x, or y) and
#' it is discrete.
ggplot(goats, aes(x=chest, y=weight, col=factor(gender))) +
  geom_point() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)") +
  scale_color_discrete("", labels=c("Male", "Female"))

#' Note that if we'd labelled things in the `data.frame` sensibly, we wouldn't
#' need to do this - let's create a new factor column with a fixed up version
#' of Gender. To do this we convert to a factor, and specify the (existing) levels
#' and what we want them labelled:
goats$Gender = factor(goats$gender, levels = c(1,2), labels=c("Male", "Female"))

#' Now the plot is simpler. In general, the advice is to keep things in the data.frame
#' using sensible labels. These may not be what ends up on a plot (as we can
#' change that if we want as above), but there is no reason to use 1 and 2 if Male
#' and Female make more sense!
ggplot(goats, aes(x=chest, y=weight, col=Gender)) +
  geom_point() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' We can change the x variable to location. Note that ggplot2 doesn't change
#' the plot type, as that is defined by the geometries we choose, and we've
#' specified `geom_point` here:
ggplot(goats, aes(x=location, y=weight, col=Gender)) +
  geom_point() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' But we can change to a boxplot easy enough:
ggplot(goats, aes(x=location, y=weight, col=Gender)) +
  geom_boxplot() +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' Note that as we've still coloured by gender, ggplot2 automatically splits all the
#' locations by gender as well. Neat!
#' 
#' An alternate would be to use **small multiples** or a faceted plot:
ggplot(goats, aes(x=location, y=weight)) +
  geom_boxplot() +
  facet_wrap(~Gender) +
  xlab("Chest girth (cm)") +
  ylab("Weight (kg)")

#' We can do facetted histograms or histograms on top of each other as well:
#' Note that `geom_histogram` uses the `x` aesthetic - the `y` value is computed
ggplot(goats, aes(x=weight, fill=Gender)) +
  geom_histogram(bins=10) +
  xlab("Weight (kg)")

#' It's easy to do side by side as well
ggplot(goats, aes(x=weight, fill=Gender)) +
  geom_histogram(bins=10) +
  facet_wrap(~Gender, ncol=1) +
  xlab("Weight (kg)")

#' Note the legend isn't much use here as we already have it split up, so
#' lets turn it off by using `guides` to set the `fill` guide to `"none"`.
#' Just for fun, I set the theme as well
ggplot(goats, aes(x=weight, fill=Gender)) +
  geom_histogram(bins=10) +
  facet_wrap(~Gender, ncol=1) +
  xlab("Weight (kg)") +
  guides(fill="none") +
  theme_bw()

#' Hopefully you see a hint of the power that ggplot2 has compared to just using
#' base graphics. That one would be quite tricky in base. I think it would be
#' done as follows by precomputing the histogram breaks?

par(mfrow=c(2,1))
male_weights = goats$weight[goats$gender == 1]
female_weights = goats$weight[goats$gender == 2]
breaks = seq(5,75,by=5) # breaks to use for histogram (common for both)
hist(male_weights, breaks=breaks, main="Male", xlab="Weight (kg)", col='red')
hist(female_weights, breaks=breaks, main="Female", xlab="Weight (kg)", col='blue')

#' but this has repetition of the x axis which isn't really needed and looks
#' a bit more naff...
#'
#' ## Exercise 2
#' 
#' ggplot2 can't do 3D (well, it sort of can, but generally there are better
#' tools such as `rgl` for that!) Nonetheless, you can do contour plots.
#' 
#' The first step here is to reshape the `volcano` data into a long data frame
#' with x, y and z. We can use `melt` out of `reshape2` library for this, or
#' probably by converting to a `data.frame`, adding a `y` column for the row
#' identifier, and then using `gather` from `tidyr`.
library(reshape2)
volcano3d <- melt(volcano, varnames=c('x', 'y'), value.name='z')
head(volcano3d)

#' Now we can use ggplot to plot it with `geom_tile`
ggplot(volcano3d, aes(x=x, y=y, fill=z)) +
  geom_tile()

#' We can add a contour to it with `geom_contour` by first adding a `z` aesthetic
ggplot(volcano3d, aes(x=x, y=y, z=z, fill=z)) +
  geom_tile() + geom_contour(col='black')

#' Colours are easy to change using `scale_fill_gradient`
#' (we want to change a scale, it's the fill aesthetic, and use
#' a gradient)
ggplot(volcano3d, aes(x=x, y=y, z=z, fill=z)) +
  geom_tile() + geom_contour(col='black') +
  scale_fill_gradient(low='white', high='red')

#' We can get rid of the extra whitespace around the figure by changing
#' the x and y scales. By default they are expanded a bit past the range
#' of the data
ggplot(volcano3d, aes(x=x, y=y, z=z, fill=z)) +
  geom_tile() + geom_contour(col='black') +
  scale_fill_gradient(low='white', high='red') +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

#' ## Exercise 3
#' 
#' Plotting linear model results for simple models can be done entirely
#' within ggplot2. You can use `geom_smooth` to add a smoother to numeric data
#' which could be a linear model (or some other type of model, such as
#' a loess smoother).
#' 
#' For data that is purely factorial like the `glaucoma` data below, you
#' need to use another trick to combine datasets in ggplot2
#' 
#' Read the data in
glaucoma = read.csv("http://www.massey.ac.nz/~jcmarsha/rcourse/glaucoma.csv")

#' Now let's visualise it
ggplot(glaucoma, aes(x=Sex, y=MD, col=Diagnosis)) +
  geom_boxplot()

#' OK, now we want to fit a linear model and then plot the results on top
#' of the boxplots. We'd first fit the linear model and create
#' a data.frame with our predictions alongside the data used for the
#' predictions:
mod = lm(MD ~ Sex*Diagnosis, data=glaucoma)
new_data = expand.grid(Sex=c("f", "m"), Diagnosis=c("G", "S"))
pred = predict(mod, new_data, interval="confidence")
pred_df = cbind(new_data, pred)
pred_df
#' OK, now we can do the plot by adding a `geom_errorbar` that uses
#' the new data source. To do this we have to move the aesthetic
#' definition into the geom rather than in the main ggplot2, as we
#' have different aesthetics for each geom. We also must use a trick
#' to get the error bars to 'dodge' correctly like the boxplot does.
#' I found the 0.75 value just by experimentation. I think you can
#' get it out of the ggplot object somehow though!
ggplot() +
  geom_boxplot(aes(x=Sex, y=MD, col=Diagnosis), data=glaucoma) +
  geom_errorbar(aes(x=Sex, ymin=lwr, ymax=upr, group=Diagnosis), position = position_dodge(width=0.75),
                width=0.2, data=pred_df) +
  geom_point(aes(x=Sex, y=fit, group=Diagnosis), position = position_dodge(width=0.75), size=3, data=pred_df)

#' Lastly, note that `visreg` can make use of `ggplot2` when it plots by adding
#' `gg=TRUE`. You can then do all the usual things with ggplot like add
#' themes etc
library(visreg)
visreg(mod, "Diagnosis", by="Sex", gg=TRUE) +
  theme_bw()

