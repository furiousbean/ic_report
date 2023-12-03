source("ic.R")

set.seed(3)
coef <- 5
s <- coef * sin(2 * pi / 6 * 1:50)
x <- s + rnorm(50)

try_ic(s, r_max = 10)
