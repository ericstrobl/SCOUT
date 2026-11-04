# SCOUT benchmark: 100 replications x 16 (n, p) settings.

library(SCOUT)
for (f in c("sim_data.R", "metrics.R", "SCOUT_no_2D.R", "SCOUT_no_COMPACT.R"))
  source(file.path("Experiments", f))

load_helper <- function(f) {
  e <- new.env(parent = globalenv())
  sys.source(file.path("Experiments", f), envir = e)
  e
}
sbw <- load_helper("sbw.R")
inv <- load_helper("invreg.R")
tmle <- load_helper("cv_tmle.R")

NSIM <- 100L
NGRID <- c(250L, 500L, 1000L, 2000L)
PGRID <- c(50L, 100L, 200L, 500L)
ALGORITHMS <- c("SCOUT", "SCOUT-no-2D", "SCOUT-no-COMPACT",
                "SBW-SCORE", "Inverse regression", "CV-TMLE")

fit_method <- function(alg, dat, seed) {
  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  X <- dat$X; T <- dat$T; Y <- dat$Y_list
  S <- visit_fits <- NULL
  target <- "X"

  if (alg %in% ALGORITHMS[1:3]) {
    # Original SCOUT CV default is seed = 1; no-COMPACT receives the benchmark seed.
    f <- switch(alg,
      "SCOUT" = SCOUT(X, T, Y),
      "SCOUT-no-2D" = SCOUT_no_2D(X, T, Y),
      "SCOUT-no-COMPACT" = SCOUT_no_COMPACT(X, T, Y, seed = seed))
    a <- f$alpha
    if (any(!is.finite(a)) || any(a < -1e-8) || abs(sum(a) - 1) > 1e-6)
      stop("Expected nonzero sum-to-one outcome weights from ", alg)
    est <- vapply(f$final, `[[`, numeric(1), "ATE")
    cate <- lapply(f$final, `[[`, "CATE")
    if (alg != "SCOUT-no-COMPACT") {
      S <- lapply(f$final, `[[`, "S_visit")
      target <- "S"
    }
  } else if (alg == "SBW-SCORE") {
    f <- sbw$sbw_score(X, T, Y, B = 0, seed = seed, verbose = FALSE)
    a <- f$alpha; est <- f$tau; cate <- f$CATE
  } else if (alg == "Inverse regression") {
    f <- inv$invreg_nn(X, T, Y, B = 0, seed = seed, verbose = FALSE)
    a <- f$alpha
    est <- vapply(Y, \(y) lm.wfit(cbind(1, T), drop(y %*% a),
                                 f$ipw)$coefficients[2], numeric(1))
    cate <- f$CATE
  } else if (alg == "CV-TMLE") {
    f <- tmle$cv_tmle_global(X, T, Y, treated = 1, higher_better = TRUE,
                            outcome_mode = "last", seed = seed, verbose = FALSE)
    a <- f$alpha_mean
    visit_fits <- lapply(seq_along(Y), \(k) tmle$tmle_fixed_alpha(
      X, T, Y[[k]], a, treated = 1, seed = seed + 1000 * k,
      return_cate = TRUE))
    est <- vapply(visit_fits, `[[`, numeric(1), "ATE")
    cate <- lapply(visit_fits, `[[`, "CATE")
  } else stop("Unknown algorithm: ", alg)

  list(alpha = a, estimate = est, CATE = cate, S = S, cate_target = target,
       runtime = proc.time()[["elapsed"]] - start, raw_fit = f,
       tmle_visit_fits = visit_fits)
}

evaluate_fit <- function(fit, dat, sim, n, p, alg) {
  K <- length(dat$Y_list)
  fit$sample_true_effect <- vapply(seq_len(K), \(k)
    mean(drop(dat$truth$tau[[k]] %*% fit$alpha)), numeric(1))
  fit$targets <- vector("list", K)
  rows <- vector("list", K)

  for (k in seq_len(K)) {
    tau_x <- drop(dat$truth$tau[[k]] %*% fit$alpha)
    tau_target <- tau_x
    retained <- NA_real_
    if (fit$cate_target == "S") {
      S <- as.matrix(fit$S[[k]])
      tau_target <- drop(lm.fit(cbind(1, S), tau_x)$fitted.values)
      retained <- var(tau_target) / var(tau_x)
    }
    fit$targets[[k]] <- list(tau_x = tau_x, tau_target = tau_target)
    err <- fit$estimate[k] - fit$sample_true_effect[k]
    rows[[k]] <- data.frame(sim = sim, n = n, p = p, algorithm = alg, visit = k,
                           bias = err, squared_error = err^2,
                           cate_squared_error = mean((fit$CATE[[k]] - tau_target)^2),
                           R2_retained = retained, runtime_sec = fit$runtime)
  }
  list(fit = fit, results = do.call(rbind, rows))
}

# Retain the exact original seed formulas and (n, p, method) iteration order.
dir.create("simulation_fits", showWarnings = FALSE)
rows <- vector("list", NSIM * length(NGRID) * length(PGRID) * length(ALGORITHMS))
row_id <- 0L
for (sim in seq_len(NSIM)) {
  message("Simulation ", sim, "/", NSIM)
  master <- make_master_data(seed = 10000 + sim)
  scenario <- 0L

  for (n in NGRID) for (p in PGRID) {
    scenario <- scenario + 1L
    dat <- get_sim_data(master, n = n, p = p, effect = "alternative")
    if (any(vapply(dat$Y_list, nrow, 0L) != nrow(dat$X)) ||
        (!is.null(dat$visit_ids) &&
         !all(vapply(dat$visit_ids, identical, logical(1), seq_len(nrow(dat$X))))))
      stop("Comparator methods require complete visits.")
    fits <- setNames(vector("list", length(ALGORITHMS)), ALGORITHMS)

    for (aa in seq_along(ALGORITHMS)) {
      alg <- ALGORITHMS[aa]
      seed <- 1000000 + 1000 * sim + 10 * scenario + aa
      fit <- fit_method(alg, dat, seed)
      fit$seed <- seed
      evaluated <- evaluate_fit(fit, dat, sim, n, p, alg)
      row_id <- row_id + 1L
      rows[[row_id]] <- evaluated$results
      fits[[alg]] <- evaluated$fit
    }
    saveRDS(list(sim = sim, n = n, p = p, data = dat, fits = fits),
            sprintf("simulation_fits/sim_%03d_n%d_p%d.rds", sim, n, p))
  }
  results <- do.call(rbind, rows[seq_len(row_id)])
  saveRDS(results, "simulation_results.rds")
  write.csv(results, "simulation_results.csv", row.names = FALSE)
}

# Visit-specific summary: genuine Monte Carlo ATE RMSE.
by <- c("n", "p", "algorithm", "visit")
summary <- aggregate(cbind(bias, squared_error, cate_squared_error, runtime_sec) ~
                     n + p + algorithm + visit, results, mean)
retained <- aggregate(R2_retained ~ n + p + algorithm + visit,
                      subset(results, !is.na(R2_retained)), mean)
summary <- merge(summary, retained, by = by, all.x = TRUE)
summary$RMSE <- sqrt(summary$squared_error)
summary$CATE_RMSE <- sqrt(summary$cate_squared_error)
summary <- summary[, c(by, "bias", "RMSE", "CATE_RMSE", "R2_retained", "runtime_sec")]
summary <- summary[order(summary$n, summary$p,
                         match(summary$algorithm, ALGORITHMS), summary$visit), ]
write.csv(summary, "simulation_summary.csv", row.names = FALSE)

results$ATE_error <- sqrt(results$squared_error)
results$CATE_RMSE <- sqrt(results$cate_squared_error)

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(Mean = NA_real_, Upper = NA_real_, Lower = NA_real_))
  m <- mean(x)
  if (n == 1L) return(c(Mean = m, Upper = NA_real_, Lower = NA_real_))
  d <- qt(0.975, df = n - 1L) * sd(x) / sqrt(n)
  c(Mean = m, Upper = m + d, Lower = m - d)
}

make_prism_table <- function(metric, n_value) {
  d <- results[results$n == n_value,
               c("sim", "n", "p", "algorithm", "visit", metric)]
  names(d)[6L] <- "value"
  per_sim <- aggregate(value ~ sim + n + p + algorithm, d,
                       \(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  out <- data.frame(X = sort(unique(per_sim$p)))
  for (alg in ALGORITHMS) {
    tmp <- per_sim[per_sim$algorithm == alg, ]
    ci <- do.call(rbind, lapply(out$X,
                               \(pp) mean_ci(tmp$value[tmp$p == pp])))
    out[[paste0(alg, " Mean")]] <- ci[, "Mean"]
    out[[paste0(alg, " Upper Limit")]] <- ci[, "Upper"]
    out[[paste0(alg, " Lower Limit")]] <- ci[, "Lower"]
  }
  out
}

# Manuscript figure metrics: ATE error, CATE RMSE, R2 retained, runtime.
dir.create("prism_tables", showWarnings = FALSE)
for (metric in c("ATE_error", "CATE_RMSE", "R2_retained", "runtime_sec"))
  for (nn in NGRID)
    write.csv(make_prism_table(metric, nn),
              sprintf("prism_tables/%s_n%d.csv", metric, nn),
              row.names = FALSE, na = "")
message("Finished; figure tables are in prism_tables/.")
