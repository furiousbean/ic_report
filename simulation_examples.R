source("simulation.R")

default_signal = 2 * sin(4 * pi / 6 * 1:50)
log_signal = 200 * log(1:50)
rank_6_signal = 2 * sin(4 * pi / 6 * 1:50) - 1.5 * sin(3 * pi / 7 * 1:50) + 1 * sin(pi / 4 * 1:50)
underest_rank4 = 2 * exp(-0.05 * 1:100) + 0.5 * exp(0.03 * 1:100) + 1.5 * cos(2 * pi * 1:100 / 30 + pi / 2)
underest2_rank7 = 0.5 * exp(0.05 * 1:100) + 7e-4 * (1:100) ^ 2 - 0.1 * 1:100 + 3 + 4.12 * cos(2 * pi * 1:100 / 4) + 4.12 * cos(2 * pi * 1:100 / 3 + pi / 2)
cores = 8
cex = 1

tryCatch({
    cores_bias <- 0
    if (.Platform$OS.type == "windows") {
        cores_bias <- 1
    }
    cores <<- max(parallel::detectCores() - cores_bias, 1)
}, error = function(e) {
    print("Please install \"parallel\" package for automatic detection of core numbers")
})

cl <- NULL

# cl <- makeCluster(getOption("cl.cores", cores))

# NKZ Examples

# pic_1_data <- different_sigmas_simul(signal = default_signal, cl = cl)
# pic_2_data <- different_sigmas_simul(signal = rank_6_signal, cl = cl)
# pic_3_data <- different_sigmas_simul(signal = log_signal, left=1e-2, right=1e4, cl = cl)
# 
# pic_4_data <- different_sigmas_simul(signal = default_signal, method = "mgn", cl = cl)
# pic_5_data <- different_sigmas_simul(signal = rank_6_signal, method = "mgn", cl = cl)
# pic_6_data <- different_sigmas_simul(signal = log_signal, left=1e-2, right=1e4, method = "mgn", cl = cl)

# save(pic_1_data, pic_2_data, pic_3_data, pic_4_data, pic_5_data, pic_6_data, file="pic.RData")

load(file="pic.RData")

plot_different_sigmas_data(pic_1_data, cex = cex)
plot_different_sigmas_data(pic_2_data, cex = cex)
plot_different_sigmas_data(pic_3_data, cex = cex)

plot_different_sigmas_data(pic_4_data, cex = cex)
plot_different_sigmas_data(pic_5_data, cex = cex)
plot_different_sigmas_data(pic_6_data, cex = cex)

plot_different_datas(pic_1_data, pic_4_data, cex = cex)
plot_different_datas(pic_2_data, pic_5_data, cex = cex)
plot_different_datas(pic_3_data, pic_6_data, cex = cex)

sample_svd <- one_sigma_simulation(default_signal, method = "svd", cl = cl)
table(sample_svd$chosen_rank_sample)
plot(sample_svd$r_net, sample_svd$reconstruction_rmse, type = "b")
sample_mgn <- one_sigma_simulation(default_signal, method = "mgn", cl = cl)
table(sample_mgn$chosen_rank_sample)
plot(sample_mgn$r_net, sample_mgn$reconstruction_rmse, type = "b")

sample_svd <- one_sigma_simulation(log_signal, method = "svd", cl = cl)
table(sample_svd$chosen_rank_sample)
plot(sample_svd$r_net, sample_svd$reconstruction_rmse, type = "b")
sample_mgn <- one_sigma_simulation(log_signal, method = "mgn", cl = cl)
table(sample_mgn$chosen_rank_sample)
plot(sample_mgn$r_net, sample_mgn$reconstruction_rmse, type = "b")

# NEG Examples

# pic_underest_data_svd <- different_sigmas_simul(signal = underest_rank4, left=1e-2, right=1e4, cl = cl)
# pic_underest_data_mgn <- different_sigmas_simul(signal = underest_rank4, left=1e-2, right=1e4, method = "mgn", cl = cl)

#save(pic_underest_data_svd, pic_underest_data_mgn, file="pic_underest.RData")

load(file="pic_underest.RData")

plot_different_sigmas_data(pic_underest_data_svd, cex = cex)
plot_different_sigmas_data(pic_underest_data_mgn, cex = cex)

plot_different_datas(pic_underest_data_svd, pic_underest_data_mgn)

# pic_underest2_data_svd <- different_sigmas_simul(signal = underest2_rank7, left=1e-2, right=1e4, cl = cl)
# pic_underest2_data_mgn <- different_sigmas_simul(signal = underest2_rank7, left=1e-2, right=1e4, method = "mgn", cl = cl)
# 
# save(pic_underest2_data_svd, pic_underest2_data_mgn, file="pic_underest2.RData")

load(file="pic_underest2.RData")

plot_different_sigmas_data(pic_underest2_data_svd, cex = cex)
plot_different_sigmas_data(pic_underest2_data_mgn, cex = cex)

plot_different_datas(pic_underest2_data_svd, pic_underest2_data_mgn, cex = cex)

# neighbors detect heuristic

sample_svd_n <- one_sigma_simulation(default_signal, method = "svd", detect = neighbors_detect, cl = cl)
hist(sample_svd_n$chosen_rank_sample)
sample_mgn_n <- one_sigma_simulation(default_signal, method = "mgn", detect = neighbors_detect, cl = cl)
hist(sample_mgn_n$chosen_rank_sample)

# sample_svd <- one_sigma_simulation(rank_6_signal, method = "svd", cl = cl, sigma = 2.069138)
# sample_svd <- one_sigma_simulation(rank_6_signal, method = "svd", cl = cl, sigma = 1.274275)
# sample_svd <- one_sigma_simulation(rank_6_signal, method = "svd", cl = cl, sigma = 0.78476)
# table(sample_svd$chosen_rank_sample)
# hist(table(sample_svd$chosen_rank_sample))

if (!is.null(cl)) {
    stopCluster(cl)
}
