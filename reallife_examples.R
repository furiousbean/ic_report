source("ic.R")

library(Rssa)
data(AustralianWine)

dat <- AustralianWine[, 1]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[, 2]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[, 3]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[, 4]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[1:174, 5]
#piece with no gaps
plot(AustralianWine[, 5])
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[, 6]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

dat <- AustralianWine[, 7]
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

#####################

data(co2)

dat <- co2
plot(dat)
res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)

r_restore <- 13
co2_signal_data <- hlra_ar(co2, r_restore)
co2_signal <- co2_signal_data$signal
co2_noise <- co2_signal_data$noise
plot(co2_signal)
plot(co2_noise)
plot(acf(co2_noise))
co2_sigma <- sd(co2_noise)

set.seed(31)
fixed_co2 <- co2_signal + rnorm(length(co2_signal)) * co2_sigma / 10
matplot(1:length(fixed_co2), cbind(co2, fixed_co2), type = "l")
res_rank_fixed <- try_ic(fixed_co2, r_max = 30, begin_idx = 2)

co2_signal_ssa <- ssa(co2_signal)
parestimate(co2_signal_ssa, groups=list(1:r_restore))

plot(ssa(fixed_co2, svd.method = "svd"))
plot(ssa(co2_signal))

#####################

library(ssabook)

data("dwarfst")

dat <- dwarfst

s <- ssa(dat)
L <- 120
plot(s)
plot(s, type = "vectors", groups = 1:30)
plot(wcor(s, groups = 1:(min(20, L))))

res_rank <- try_ic(dat, r_max = 30, begin_idx = 3)

#####################

data("USUnemployment")

dat <- USUnemployment[, 3] # works on it
# dat <- USUnemployment[, 2] # heavy red noise -- it dies
s <- ssa(dat)
L <- 204
plot(s)
plot(s, type = "vectors", groups = 1:30)
plot(wcor(s, groups = 1:(min(20, L))))

res_rank <- try_ic(dat, r_max = 30, begin_idx = 2)
