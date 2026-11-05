
# SCOUT CRT simulations: three conditions, three visits, mean ATE / mean Q
source("sim_data.R"); source("randomization_dgp.R"); source("crt_scout_fast.R")
M <- 500L; B <- 500L; n <- 1000L; p <- 200L
DGP_SEED <- sample.int(.Machine$integer.max, 1L); VERSION <- 1L
CRT_CORES <- 4L; CRT_RESELECT_ETA <- TRUE
LEVEL <- .05; MIN_VALID <- .90; PROP_CLIP <- 1e-6
OUT <- sprintf("scout_CRT_M%d_B%d_n%d_p%d_seed%d_v%d", M, B, n, p, DGP_SEED, VERSION)
dir.create(OUT, showWarnings = FALSE)
cat("DGP seed:", DGP_SEED, "\n")

Q_stat <- function(z) sqrt(mean((z - mean(z))^4))
p_upper <- function(x, z) (1 + sum(z >= x))/(length(z) + 1)
extract_stats <- function(f) {
  a <- vapply(f$final, `[[`, 0.0, "ATE")
  C <- do.call(cbind, lapply(f$final, `[[`, "CATE"))
  q <- apply(C, 2, Q_stat)
  list(ATE = c(a, mean(a)), Q = c(q, mean(q)))
}
make_constant <- function(d, scales, delta) {
  for (k in seq_along(d$Y_list))
    d$Y_list[[k]] <- d$Y_list[[k]] + tcrossprod(d$T, scales[k] * delta)
  d
}

crt_one <- function(b, s) {
  set.seed(s$seed + 100000L * (s$K + 1L) + b)
  ts <- rbinom(s$n, 1L, s$e)
  empty <- list(ATE = rep(NA_real_, s$J), Q = rep(NA_real_, s$J))
  if (min(tabulate(ts + 1L, nbins = 2L)) < 5L) return(empty)
  one <- function(Y) tryCatch(crt_scout_multi(s$X, ts, list(Y), eta = s$eta,
                                              seed = 1)[[1]], error = \(e) NULL)
  stats <- function(f) if (is.null(f)) empty else
    tryCatch(extract_stats(f), error = \(e) empty)
  if (s$zero) return(stats(one(s$Y)))
  Yh <- lapply(seq_len(s$K), \(k)
    s$Y[[k]] + tcrossprod(ts - s$T, s$delta[k, ]))
  pair <- tryCatch(crt_scout_multi(s$X, ts, list(s$Y, Yh), eta = s$eta,
                                   seed = 1), error = \(e) NULL)
  if (is.null(pair)) pair <- list(one(s$Y), one(Yh))
  list(ATE = stats(pair[[1]])$ATE, Q = stats(pair[[2]])$Q)
}
crt_worker <- function(b) crt_one(b, .crt_state)

run_condition <- function(d, seed, zero, cl) {
  X <- as.matrix(d$X); T <- as.numeric(d$T)
  Y <- lapply(d$Y_list, as.matrix); K <- length(Y); J <- K + 1L
  visits <- c(paste0("visit", seq_len(K)), "mean_visit_Q")
  fit <- crt_scout(X, T, Y, seed = 1, estimate_item_effects = !zero)
  obs <- extract_stats(fit)
  S <- do.call(cbind, lapply(fit$final, `[[`, "S"))
  e <- glm.fit(cbind(1, S), T, family = binomial())$fitted.values
  e <- pmin(pmax(e, PROP_CLIP), 1 - PROP_CLIP)
  state <- list(X = X, T = T, Y = Y, e = e, n = nrow(X), K = K, J = J,
                seed = seed, zero = zero, delta = fit$constant_item_effects,
                eta = if (CRT_RESELECT_ETA) NULL else fit$eta)
  if (is.null(cl)) {
    draws <- lapply(seq_len(B), crt_one, s = state)
  } else {
    parallel::clusterCall(cl, function(s) {
      assign(".crt_state", s, envir = .GlobalEnv); NULL
    }, state)
    draws <- parallel::parLapplyLB(cl, seq_len(B), crt_worker, chunk.size = 5L)
  }
  A <- do.call(rbind, lapply(draws, `[[`, "ATE"))
  Q <- do.call(rbind, lapply(draws, `[[`, "Q"))
  colnames(A) <- colnames(Q) <- visits
  rows <- lapply(seq_len(J), \(j) {
    a <- A[is.finite(A[, j]), j]; q <- Q[is.finite(Q[, j]), j]
    if (min(length(a), length(q)) < MIN_VALID * B)
      stop("Too few valid CRT draws for ", visits[j])
    pa <- p_upper(obs$ATE[j], a); pq <- p_upper(obs$Q[j], q)
    data.frame(visit = visits[j], ATE = obs$ATE[j], CATE_Q = obs$Q[j],
               ATE_p_positive = pa, CATE_p_homogeneity = pq,
               ATE_reject = pa <= LEVEL, CATE_reject = pq <= LEVEL)
  })
  list(results = do.call(rbind, rows),
       raw = list(ATE_null = A, Q_homogeneity_null = Q))
}

# Single persistent PSOCK cluster; share data once per condition.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
cl <- if (CRT_CORES > 1L) parallel::makeCluster(CRT_CORES) else NULL
if (!is.null(cl)) {
  parallel::clusterExport(cl, c("crt_scout_multi", "Q_stat", "extract_stats",
                                "crt_one", "crt_worker"), envir = .GlobalEnv)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE))
      RhpcBLASctl::blas_set_num_threads(1L)
    NULL
  })
}

dgp <- make_dgp(DGP_SEED)
saveRDS(dgp, file.path(OUT, "fixed_dgp.rds"))
all_results <- vector("list", M)
tryCatch({
  for (sim in seq_len(M)) {
    message("Replication ", sim, "/", M)
    path <- file.path(OUT, sprintf("rep_%03d.rds", sim))
    if (file.exists(path)) {
      obj <- readRDS(path)
    } else {
      master <- sample_from_dgp(dgp, seed = 10001L + sim, n_max = n)
      d0 <- get_sim_data(master, n, p, "null")
      d1 <- get_sim_data(master, n, p, "alternative")
      set.seed(300001L + sim)
      delta <- runif(ncol(d0$Y_list[[1L]]), -1, 1)
      cases <- list(no_effect = d0,
                    constant_effect = make_constant(d0, d1$truth$visit_scale, delta),
                    heterogeneous_effect = d1)
      seed <- 200001L + (B + 1L) * sim
      obj <- lapply(names(cases), \(e)
        run_condition(cases[[e]], seed, zero = (e == "no_effect"), cl = cl))
      names(obj) <- names(cases)
      saveRDS(obj, path)
    }
    all_results[[sim]] <- do.call(rbind, lapply(names(obj), \(e) {
      z <- obj[[e]]$results; z$effect <- e; z$sim <- sim; z
    }))
    if (sim %% 10L == 0L || sim == M)
      write.csv(do.call(rbind, all_results[seq_len(sim)]),
                file.path(OUT, "results_updated.csv"), row.names = FALSE)
  }
}, finally = { if (!is.null(cl)) parallel::stopCluster(cl) })
results <- do.call(rbind, all_results)
saveRDS(results, file.path(OUT, "results_updated.rds"))

# Prism-ready calibration curves, Type I error, and power (exact 95% CIs).
visits <- unique(results$visit); grid <- seq(.01, 1, .01)
curve <- function(d, col) 100 * vapply(grid, \(a) mean(d[[col]] <= a), 0.0)
for (v in visits) {
  no <- subset(results, effect == "no_effect" & visit == v)
  co <- subset(results, effect == "constant_effect" & visit == v)
  panel <- data.frame(Nominal_Percent = 100 * grid,
    ATE_ZeroEffect_Percent = curve(no, "ATE_p_positive"),
    CATE_ConstantEffect_Percent = curve(co, "CATE_p_homogeneity"))
  fname <- if (v == "mean_visit_Q") "Prism_Panel_A_Calibration_updated.csv" else
    paste0("Prism_Panel_A_Calibration_", v, ".csv")
  write.csv(panel, file.path(OUT, fname), row.names = FALSE)
}
rate_ci <- function(x) c(mean(x), unname(binom.test(sum(x), length(x))$conf.int))
spec <- function(e, v, col, name) data.frame(e, v, col, name)
report <- function(s) do.call(rbind, lapply(seq_len(nrow(s)), \(i) {
  d <- subset(results, effect == s$e[i] & visit == s$v[i])
  z <- 100 * rate_ci(d[[s$col[i]]] <= LEVEL)
  data.frame(Test = s$name[i], Rate = z[1], Lower = z[2], Upper = z[3],
             Error_Minus = z[1] - z[2], Error_Plus = z[3] - z[1])
}))
TypeI <- do.call(rbind, lapply(visits, \(v) rbind(
  spec("no_effect", v, "ATE_p_positive", paste("ATE null -", v)),
  spec("constant_effect", v, "CATE_p_homogeneity", paste("CATE constant -", v)))))
Power <- do.call(rbind, lapply(visits, \(v) rbind(
  spec("constant_effect", v, "ATE_p_positive", paste("ATE constant -", v)),
  spec("heterogeneous_effect", v, "ATE_p_positive", paste("ATE heterogeneous -", v)),
  spec("heterogeneous_effect", v, "CATE_p_homogeneity", paste("CATE heterogeneous -", v)))))
write.csv(report(TypeI), file.path(OUT, "Prism_Panel_B_TypeI_updated.csv"), row.names = FALSE)
write.csv(report(Power), file.path(OUT, "Prism_Panel_C_Power_updated.csv"), row.names = FALSE)
cat("Complete:", M, "replications x", B, "randomizations. Output:", OUT, "\n")
