#compute the power for a clustered design
powa <- function(n, n2=n, p=0.45, delta_p, deff=2, alpha=0.05) {
  p2 <- p + delta_p
  q1 <- 1 - p
  q2 <- 1 - p2
  pm <- (n * p + n2 * p2)/(n + n2)
  z <- qnorm(1 - alpha/2)
  ds <- z * sqrt(deff*(1/n + 1/n2) * pm * (1 - pm))
  ex <- abs(delta_p)
  sd <- sqrt(deff*(p * q1/n + p2 * q2/n2))
  1 - pnorm((ds - ex)/sd) + pnorm((-ds - ex)/sd)
}

# this will be zero when the power given by the design matches the power we pass in
powa_match <- function(delta_p, power, n=300, p=0.45, deff=2, alpha=0.05) {
  powa(n=n, p=p, delta_p=delta_p, deff=deff, alpha=alpha) - power
}

# now find the delta_p that gives the required power.
# idea is we need to find the delta_p such that powa_match is zero. The `uniroot`\
# function in R does that (finds a 'root' of a function)
find_delta <- function(power, n=300, p=0.45, deff=2, alpha=0.05) {
  uniroot(powa_match, interval=c(1e-3, 0.5), n=n, p=p, power=power, deff=deff, alpha=alpha)$root
}

# e.g. what can we detect for a power of 80% in a design with n=300, p=0.3, deff=2?
find_delta(power=0.8, n=300, p=0.3, deff=1.5)
# seems we can detect a difference in prevalence of 0.156

# with this we can generate things really fast
pow_range <- expand.grid(n=seq(200,400,by=100), p=seq(0.01,0.4,by=0.01), power=c(0.6,0.8),
                         Deff=seq(1.2,2.4,by=0.4), alpha=c(0.01, 0.05))

# for each row, compute the power
for (i in 1:nrow(pow_range)) {
  p <- pow_range[i,]
  pow_range$delta_p[i] = find_delta(power=p$power, n=p$n, p=p$p, deff=p$Deff, alpha=p$alpha)
}

# and we can do a plot
library(ggplot2)
ggplot(pow_range, aes(x=p, y=delta_p, col=factor(power), lty=factor(alpha))) +
  geom_line() +
  facet_grid(Deff~n)
