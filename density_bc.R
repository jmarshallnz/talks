# code for a boundary corrected density estimate
density_bc <- function(x, x0 = -Inf, x1 = Inf)
{
  d <- density(x)
  inside <- d$x > x0 && d$x < x0
  d$y <- d$y / pnorm(x1, d$x, d$bw) - pnorm(x0, d$x, d$bw)
  d$x <- d$x[inside]
  d$y <- d$y[inside]
  d
}
