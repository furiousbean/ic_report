source("ic.R")
library(snow)

default_t = 100
colors <- c("black", "red", "darkgreen", "darkblue", "violet", "yellow3")

failed_series <- NULL
failed_rank <- NULL

primitive_detect <- function(ic_data, critetion = "bic") {
    index <- as.numeric(which.max(ic_data[, critetion]))
    ic_data[index, "r"][[1]]
}

neighbors_detect <- function(ic_data, critetion = "bic") {
    critretion_data <- as.numeric(ic_data[, critetion])
    if (length(critretion_data) == 1) {
        return(ic_data[1, "r"][[1]])
    }

    if (length(critretion_data) == 2) {
        if (critretion_data[1] >= critretion_data[2]) {
            return(ic_data[1, "r"][[1]])
        } else {
            return(ic_data[2, "r"][[1]])
        }
    }

    for (i in 2:(length(critretion_data) - 1)) {
        if (critretion_data[i] >= critretion_data[i - 1] &&
            critretion_data[i] >= critretion_data[i + 1]) {
            return(ic_data[i, "r"][[1]])
        }
    }

    m <- length(critretion_data)
    if (critretion_data[1] >= critretion_data[m]) {
        return(ic_data[1, "r"][[1]])
    } else {
        return(ic_data[m, "r"][[1]])
    }
}

one_sigma_simulation <- function(signal, t = default_t, method = "svd", criterion = "bic",
    detect = primitive_detect, seed = 15,
    sigma = 1,
    simul = function(...) signal + sigma * rnorm(length(signal)),
    r_range = 0:10,
    L = NULL,
    cl = NULL) {
    set.seed(seed)
    data <- lapply(1:t, simul)

    lapply_fun <- lapply

    if (!is.null(cl)) {
        lapply_fun <- function(...) parLapply(cl, ...)
        clusterCall(cl, function() {
            library(svd)
            library(Matrix)
            library(rhlra)
            library(Rssa)
        })
        clusterExport(cl, c("calculate_ic_factory",
                            "primitive_detect",
                            "neighbors_detect",
                            "eval_ic"))
    }

    all_data <- lapply_fun(1:t, function(i) {
        cat(method, " ", i, "\n")
        series <- data[[i]]
        f <- calculate_ic_factory(method = method, L = L, r_range = r_range)
        ic_data <- f(series, signal)
        best_rank <- detect(ic_data, criterion)
        r_net <- as.numeric(ic_data$r)
        signal_dist_ic_hlra <- NA
        if (method == "svd") {
          call_pars <- list(series = series, r = best_rank)
          if (!is.null(L)) {
            call_pars$L <- L
          }
          tryCatch(
            expr = {
              if (best_rank == 0) {
                signal_dist_ic_hlra <-  mean(signal^2)
              } else {
                hlra_obj <- do.call("hlra", call_pars)
                signal_dist_ic_hlra <-  mean((signal - hlra_obj$signal)^2)
              }
            },
            error = function(e){ 
              failed_series <<- series
              failed_rank <<- best_rank
              failed_L <<- L
              print(r_net)
              print(best_rank)
            }
          )
        }
        answer <- list(best_rank = best_rank,
                       signal_dist = as.numeric(ic_data$signal_dist),
                       signal_dist_ic = as.numeric(ic_data$signal_dist)[r_net == best_rank],
                       signal_dist_ic_hlra = signal_dist_ic_hlra,
                       full_dist = mean((series - signal)^2))
        if (i == 1) {
            answer$r_net <- r_net
        }
        answer
    })

    r_net <- all_data[[1]]$r_net
    rank_data <- sapply(all_data, function(x) x$best_rank)
    signal_dist_sums <- rowSums(sapply(all_data, function(x) x$signal_dist))
    signal_dist_ic_sum <- sum(sapply(all_data, function(x) x$signal_dist_ic))
    signal_dist_ic_hlra_sum <- sum(sapply(all_data, function(x) x$signal_dist_ic_hlra))
    full_dist_sum <- sum(sapply(all_data, function(x) x$full_dist))

    best_reconstruction_rank <- r_net[which.min(signal_dist_sums)]
    list(chosen_rank_sample = rank_data,
         best_reconstruction_rank = best_reconstruction_rank,
         r_net = r_net,
         reconstruction_rmse = sqrt(signal_dist_sums / t),
         best_reconstruction_rmse = sqrt(min(signal_dist_sums) / t),
         chosen_rank_rmse = sqrt(signal_dist_ic_sum / t),
         full_dist_rmse = sqrt(full_dist_sum / t),
         hlra_chosen_rank_rmse = sqrt(signal_dist_ic_hlra_sum / t)
         )
}

different_sigmas_simul <- function(signal = default_signal,
                                   left = 1e-2, right = 1e2, sigmas_length = 20,
                                   sigmas = exp(seq(log(left), log(right), length.out = sigmas_length)),
                                   method = "svd", criterion = "bic", L = NULL, cl = NULL) {
    best_r_rank <- numeric(0)
    avg_chosen_rank <- numeric(0)
    best_r_rmse <- numeric(0)
    chosen_r_rmse <- numeric(0)
    hlra_chosen_r_rmse <- numeric(0)
    full_dist_rmse <- numeric(0)
    sapply(sigmas, function(sigma) {
        cat("sigma = ", sigma, "\n")
        info <- one_sigma_simulation(signal = signal, sigma = sigma, method = method, criterion = criterion,
                                     L = L, cl = cl)
        best_r_rank <<- c(best_r_rank, info$best_reconstruction_rank)
        avg_chosen_rank <<- c(avg_chosen_rank, mean(info$chosen_rank_sample))
        best_r_rmse <<- c(best_r_rmse, info$best_reconstruction_rmse)
        chosen_r_rmse <<- c(chosen_r_rmse, info$chosen_rank_rmse)
        full_dist_rmse <<- c(full_dist_rmse, info$full_dist_rmse)
        hlra_chosen_r_rmse <<- c(hlra_chosen_r_rmse, info$hlra_chosen_rank_rmse)
    })
    list(sigmas = sigmas, best_r_rmse = best_r_rmse, chosen_r_rmse = chosen_r_rmse, full_dist_rmse = full_dist_rmse,
         best_r_rank = best_r_rank, avg_chosen_rank = avg_chosen_rank, hlra_chosen_r_rmse = hlra_chosen_r_rmse)
}

plot_different_sigmas_data <- function(data, mode = NULL, cex = 1) {
    if (is.null(mode) || mode == 1) {
        matplot(data$sigmas, cbind(data$best_r_rmse, data$chosen_r_rmse),
                log="x", type = "b", pch = 1:2, col = colors[1:2], xlab = expression(sigma), ylab = "RMSE",
                cex = cex)
        legend("topleft", c("r с наименьшим RMSE", "r, выбранный по IC"), pch = 1:2, col = colors[1:2], cex = cex)
    }

    if (is.null(mode) || mode == 2) {
        matplot(data$sigmas,
                cbind(data$best_r_rmse/data$sigmas, data$chosen_r_rmse/data$sigmas, data$full_dist_rmse/data$sigmas),
                log="xy", type = "b", pch = 1:3, col = colors[1:3], xlab = expression(sigma),
                ylab = paste("RMSE относительно", expression(sigma)), cex = cex)
        legend("bottomleft", c("r с наименьшим RMSE", "r, выбранный по IC", "Максимальный r"),
               pch = 1:3, col = colors[1:3], cex = cex)
    }

    if (is.null(mode) || mode == 3) {
        matplot(data$sigmas, cbind(data$best_r_rank, data$avg_chosen_rank), log="x", type = "b",
                pch = 1:2, col = colors[1:2], xlab = expression(sigma), ylab = "r", cex = cex)
        legend("topright", c("r с наименьшим RMSE", "Средний r, выбранный по IC"), pch = 1:2, col = colors[1:2], cex = cex)
    }
}

plot_different_datas <- function(data_svd, data_mgn, cex = 1) {
    matplot(data_svd$sigmas,
            cbind(data_svd$best_r_rmse/data_svd$sigmas,
                  data_svd$chosen_r_rmse/data_svd$sigmas,
                  data_mgn$best_r_rmse/data_svd$sigmas,
                  data_mgn$chosen_r_rmse/data_svd$sigmas,
                  data_svd$hlra_chosen_r_rmse/data_svd$sigmas,
                  data_svd$full_dist_rmse/data_svd$sigmas),
            log="xy", type = "b", pch = 1:6, col = colors[1:6], xlab = expression(sigma),
            ylab = paste("RMSE относительно", expression(sigma)), cex = cex)
    legend("bottomleft", c("r с наименьшим RMSE, оценка SSA", "r, выбранный по IC(SVD), оценка SSA",
                           "r с наименьшим RMSE, оценка MGN", "r, выбранный по IC(MGN), оценка MGN",
                           "r, выбранный по IC(SVD), оценка MGN", "Максимальный r (оценка = зашумл. ряд)"),
           pch = 1:6, col = colors[1:6], cex = cex)
}
