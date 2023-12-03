library(Rssa)
library(rhlra)

eval_ic <- function(noise, r, df, n_mean = length(noise), n = length(noise)) {
    noise_norm_sq <- sum(noise^2)
    loglikelihood <- -n / 2 * log(noise_norm_sq / n_mean)
    bic <- -log(n) * df  + 2 * loglikelihood
    aic <- -2 * df  + 2 * loglikelihood
    list(r = r, df = df, loglikelihood = loglikelihood, bic = bic, aic = aic)
}

calculate_ic_factory <- function(method = "ssa", L = NULL, r_range = 0:10) {
    if (method == "ssa") {
        stop("Unjustified method")
        # return(function(series) {
        #     call_pars <- list(x = series)
        #     if (!is.null(L)) {
        #         call_pars$L <- L
        #     }
        #     ssaobj <- do.call("ssa", call_pars)
        #     ics <- sapply(r_range, function(r) {
        #         if (r == 0) {
        #             return(eval_ic(series, 0, 0))
        #         }
        #
        #         reconstrunction <- reconstruct(ssaobj, list(1:r))[[1]]
        #         noise <- series - reconstrunction
        #         # eval_ic(noise, r, 2 * r)
        #         eval_ic(noise, r, 26 / 9 * r) # подгон
        #     })
        #     data.frame(t(ics))
        # });
    }

    if (method == "svd") {
        return(function(series, signal = NULL) {
            call_pars <- list(x = series, svd.method = "svd")
            # SVD is slow but accurate
            if (!is.null(L)) {
                call_pars$L <- L
            }
            ssaobj <- do.call("ssa", call_pars)
            sigmas <- ssaobj$sigma
            L_real <- length(ssaobj$U[, 1])
            K <- length(series) - L_real + 1

            ics <- sapply(r_range, function(r) {
                if (r == 0) {
                    answer <- eval_ic(sigmas, 0, 0, n_mean = L_real * K,
                                      n = length(series))
                    if (!is.null(signal)) {
                        answer$signal_dist <- mean(signal^2)
                    }
                    return(answer)
                }

                pseudo_noise <- sigmas
                pseudo_noise[1:r] <- 0

                answer <- eval_ic(pseudo_noise, r,
                                  length(series) * (r * (L_real + K) - r^2) / (L_real * K),
                                  n_mean = L_real * K, n = length(series))
                if (!is.null(signal)) {
                    answer$signal_dist <- mean((signal - reconstruct(ssaobj, list(1:r))[[1]])^2)
                }
                answer
            })
            data.frame(t(ics))
        });
    }

    if (method == "mgn") {
        return(function(series, signal = NULL) {
            ics <- sapply(r_range, function(r) {
                if (r == 0) {
                    answer <- eval_ic(series, 0, 0)
                    if (!is.null(signal)) {
                        answer$signal_dist <- mean(signal^2)
                    }
                    return(answer)
                }

                call_pars <- list(series = series, r = r)
                if (!is.null(L)) {
                    call_pars$L <- L
                }
                hlra_obj <- do.call("hlra", call_pars)
                answer <- eval_ic(hlra_obj$noise, r, 2 * r)
                if (!is.null(signal)) {
                    answer$signal_dist <- mean((signal - hlra_obj$signal)^2)
                }
                answer
            })
            data.frame(t(ics))
        });
    }

    stop("Unknown method")
}

try_ic <- function(x, L = NULL, r_max = 20, begin_idx = 0, method = "svd") {
  s <- x
  r_range <- 0:r_max
  h <- calculate_ic_factory(method = method, r_range = r_range, L = L)
  svd_result <- h(x, s)
  print(r_range[which.max(svd_result$bic)])

  matplot(svd_result$r[begin_idx:length(svd_result$r)],
          cbind(svd_result$aic, svd_result$bic)[begin_idx:length(svd_result$r),],
          pch = 1:2, xlab = "rank", ylab = "IC value")
  legend("bottomright", pch = 1:2, legend = c("AIC", "BIC"), col = c(1:2))
  r_range[which.max(svd_result$bic)]
}
