# SCOUT simulation benchmark: 100 replications x 16 (n, p) settings.

SRC_DIR <- "R"  # Change to "." if the scripts are in the working directory.
for (file in c("sim_data.R", "metrics.R", "SCOUT.R",
               "SCOUT_no_2D.R", "SCOUT_no_COMPACT.R")) {
  source(file.path(SRC_DIR, file))
}

# Isolate comparator helpers that define functions with overlapping names.
load_helper <- function(file) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(SRC_DIR, file), envir = env)
  env
}
sbw <- load_helper("sbw.R")
inv <- load_helper("invreg.R")
tmle <- load_helper("cv_tmle.R")

NSIM <- 100L
NGRID <- c(250L, 500L, 1000L, 2000L)
PGRID <- c(50L, 100L, 200L, 500L)
ALGORITHMS <- c("SCOUT", "SCOUT-no-2D", "SCOUT-no-COMPACT",
                "SBW-SCORE", "Inverse regression", "CV-TMLE")

# Fit the complete native method, preserving its learned outcome and CATE.
fit_method <- function(alg, dat, seed) {
  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  X <- dat$X; T <- dat$T; Y <- dat$Y_list
  S <- visit_fits <- NULL
  target <- "X"

  if (alg %in% ALGORITHMS[1:3]) {
    # All three functions already normalize their NNLS weights.
    f <- switch(alg,
      "SCOUT" = SCOUT(X, T, Y),
      "SCOUT-no-2D" = SCOUT_no_2D(X, T, Y),
      "SCOUT-no-COMPACT" = SCOUT_no_COMPACT(X, T, Y, seed = seed)
    )
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
    # Fixed learned outcome; weighted ATE at each visit.
    est <- vapply(Y, \(y) lm.wfit(cbind(1, T), drop(y %*% a),
                                 f$ipw)$coefficients[2], numeric(1))
    cate <- f$CATE
  } else if (alg == "CV-TMLE") {
    f <- tmle$cv_tmle_global(X, T, Y, treated = 1, higher_better = TRUE,
                            outcome_mode = "last", seed = seed, verbose = FALSE)
    a <- f$alpha_mean
    # Apply the same learned outcome to each visit's TMLE.
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

# Score-based methods are evaluated against the projected treatment effect;
# other methods target the full X-indexed treatment effect.
evaluate_fit <- function(fit, dat, sim, n, p, alg) {
  K <- length(dat$Y_list)
  truth <- vapply(seq_len(K), \(k) true_effect(fit$alpha, dat$truth,
                                               visit = k), numeric(1))
  sample_truth <- vapply(seq_len(K), \(k)
    mean(drop(dat$truth$tau[[k]] %*% fit$alpha)), numeric(1))
  fit$true_effect <- truth
  fit$sample_true_effect <- sample_truth
  fit$targets <- vector("list", K)

  rows <- vector("list", K)
  for (k in seq_len(K)) {
    tau_x <- drop(dat$truth$tau[[k]] %*% fit$alpha)
    tau_target <- tau_x
    retained <- dim_score <- NA_real_
    if (fit$cate_target == "S") {
      S <- as.matrix(fit$S[[k]])
      tau_target <- drop(lm.fit(cbind(1, S), tau_x)$fitted.values)
      retained <- var(tau_target) / var(tau_x)
      dim_score <- ncol(S)
    }
    fit$targets[[k]] <- list(tau_x = tau_x, tau_target = tau_target)
    cate <- fit$CATE[[k]]
    cate_cor <- if (sd(cate) > .Machine$double.eps &&
                    sd(tau_target) > .Machine$double.eps)
      cor(cate, tau_target) else NA_real_
    bias <- fit$estimate[k] - sample_truth[k]

    rows[[k]] <- data.frame(sim = sim, n = n, p = p, algorithm = alg, visit = k,
               true_effect = truth[k], sample_true_effect = sample_truth[k],
               estimate = fit$estimate[k], bias = bias, squared_error = bias^2,
               cate_squared_error = mean((cate - tau_target)^2),
               CATE_cor = cate_cor, R2_retained = retained,
               score_dim = dim_score, runtime_sec = fit$runtime)
  }
  list(fit = fit, results = do.call(rbind, rows))
}

# Run simulations
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
      stop("This benchmark requires complete visits; the comparator methods ",
           "do not yet support visit_ids.")
    saved_fits <- setNames(vector("list", length(ALGORITHMS)), ALGORITHMS)

    for (aa in seq_along(ALGORITHMS)) {
      alg <- ALGORITHMS[aa]
      seed <- 1000000 + 1000 * sim + 10 * scenario + aa
      fit <- fit_method(alg, dat, seed)
      fit$seed <- seed
      evaluated <- evaluate_fit(fit, dat, sim, n, p, alg)
      row_id <- row_id + 1L
      rows[[row_id]] <- evaluated$results
      saved_fits[[alg]] <- evaluated$fit
    }

    # Keep fitted models and data for reproducibility and post-hoc analyses.
    saveRDS(list(sim = sim, n = n, p = p, data = dat, fits = saved_fits),
            sprintf("simulation_fits/sim_%03d_n%d_p%d.rds", sim, n, p))
  }

  # Checkpoint after each replication.
  results <- do.call(rbind, rows[seq_len(row_id)])
  saveRDS(results, "simulation_results_normalized.rds")
  write.csv(results, "simulation_results_normalized.csv", row.names = FALSE)
}

# Main summary
by <- c("n", "p", "algorithm", "visit")
summary <- aggregate(cbind(true_effect, bias, squared_error, runtime_sec) ~
                       n + p + algorithm + visit, results, mean)
cate_mse <- aggregate(cate_squared_error ~ n + p + algorithm + visit,
                      results, mean)
cate_cor <- aggregate(CATE_cor ~ n + p + algorithm + visit,
                      subset(results, !is.na(CATE_cor)), mean)
retained <- aggregate(cbind(R2_retained, score_dim) ~ n + p + algorithm + visit,
                      subset(results, !is.na(R2_retained)), mean)
for (part in list(cate_mse, cate_cor, retained))
  summary <- merge(summary, part, by = by, all.x = TRUE)
summary$RMSE <- sqrt(summary$squared_error)
summary$CATE_RMSE <- sqrt(summary$cate_squared_error)
summary$squared_error <- summary$cate_squared_error <- NULL
summary <- summary[, c(by, "true_effect", "bias", "RMSE", "CATE_RMSE",
                       "CATE_cor", "R2_retained", "score_dim", "runtime_sec")]
summary <- summary[order(summary$n, summary$p,
                         match(summary$algorithm, ALGORITHMS), summary$visit), ]
write.csv(summary, "simulation_summary.csv", row.names = FALSE)
