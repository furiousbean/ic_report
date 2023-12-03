source("ic.R")
source("simulation.R")

cex <- 1
cores = 8


set.seed(3)
coef <- 5
s <- coef * sin(2 * pi / 6 * 1:100)
x <- s + rnorm(100)
s[10:20] <- NA


try_ic(s, r_max = 10)

tryCatch({
    cores_bias <- 0
    if (.Platform$OS.type == "windows") {
        cores_bias <- 1
    }
    cores <<- max(parallel::detectCores() - cores_bias, 1)
}, error = function(e) {
    print("Please install \"parallel\" package for automatic detection of core numbers")
})

# cl <- NULL

cl <- makeCluster(getOption("cl.cores", cores))

pic_data <- different_sigmas_simul(signal = s, cl = cl, left=1e-1)

plot_different_sigmas_data(pic_data, cex = cex)

if (!is.null(cl)) {
    stopCluster(cl)
}
