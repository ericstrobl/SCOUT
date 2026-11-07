
# SCOUT CRT simulations: three conditions, three visits, mean ATE / mean Q
source("sim_data.R")
source("randomization_dgp.R")
source("crt_scout_fast.R")

M <- 500L; B <- 500L; n <- 1000L; p <- 200L
DGP_SEED <- 123456789L
CRT_CORES <- 4L
LEVEL <- .05; MIN_VALID <- .90; PROP_CLIP <- 1e-6
OUT <- sprintf("scout_CRT_M%d_B%d_n%d_p%d_dgp%d", M, B, n, p, DGP_SEED)
dir.create(OUT, showWarnings = FALSE)

Q_stat <- function(z) sqrt(mean((z - mean(z))^4))
p_upper <- function(x, z) (1 + sum(z >= x)) / (length(z) + 1)

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

  one <- function(Y) tryCatch(
    crt_scout_multi(s$X, ts, list(Y), eta = NULL, seed = 1L)[[1L]],
    error = function(e) NULL)
  stats <- function(f) if (is.null(f)) empty else
    tryCatch(extract_stats(f), error = function(e) empty)

  if (s$zero) return(stats(one(s$Y)))

  Yh <- lapply(seq_len(s$K), function(k)
    s$Y[[k]] + tcrossprod(ts - s$T, s$delta[k, ]))
  pair <- tryCatch(crt_scout_multi(s$X, ts, list(s$Y, Yh),
                                    eta = NULL, seed = 1L),
                   error = function(e) NULL)
  if (is.null(pair)) pair <- list(one(s$Y), one(Yh))
  list(ATE = stats(pair[[1L]])$ATE, Q = stats(pair[[2L]])$Q)
}

crt_worker <- function(b) crt_one(b, .crt_state)

run_condition <- function(d, seed, zero, cl) {
  X <- as.matrix(d$X); T <- as.numeric(d$T)
  Y <- lapply(d$Y_list, as.matrix); K <- length(Y); J <- K + 1L
  stopifnot(K == 3L)
  visits <- c(paste0("visit", seq_len(K)), "mean")

  fit <- crt_scout(X, T, Y, seed = 1L, estimate_item_effects = !zero)
  obs <- extract_stats(fit)
  if (any(!is.finite(unlist(obs))))
    stop("Nonfinite observed SCOUT statistics")

  S <- do.call(cbind, lapply(fit$final, `[[`, "S"))
  e <- suppressWarnings(glm.fit(cbind(1, S), T,
                                family = binomial())$fitted.values)
  if (any(!is.finite(e))) stop("Nonfinite fitted treatment propensity")
  e <- pmin(pmax(e, PROP_CLIP), 1 - PROP_CLIP)

  delta <- fit$constant_item_effects
  if (!zero && (!identical(dim(delta), c(K, ncol(Y[[1L]]))) ||
                any(!is.finite(delta))))
    stop("Invalid estimated constant outcome effects")

  state <- list(X = X, T = T, Y = Y, e = e, n = nrow(X), K = K, J = J,
                seed = seed, zero = zero, delta = delta)

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

  do.call(rbind, lapply(seq_len(J), function(j) {
    a <- A[is.finite(A[, j]), j]
    q <- Q[is.finite(Q[, j]), j]
    if (min(length(a), length(q)) < MIN_VALID * B)
      stop("Too few valid CRT draws for ", visits[j])
    data.frame(visit = visits[j], ATE = obs$ATE[j], CATE_Q = obs$Q[j],
               ATE_p_positive = p_upper(obs$ATE[j], a),
               CATE_p_homogeneity = p_upper(obs$Q[j], q))
  }))
}

# Persistent PSOCK cluster; share data once per condition.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
cl <- if (CRT_CORES > 1L) parallel::makeCluster(CRT_CORES) else NULL

dgp <- make_dgp(DGP_SEED)
saveRDS(dgp, file.path(OUT, "fixed_dgp.rds"))
all_results <- vector("list", M)

tryCatch({
  if (!is.null(cl)) {
    parallel::clusterExport(cl, c("crt_scout_multi", "Q_stat",
                                  "extract_stats", "crt_one", "crt_worker"),
                            envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_NUM_THREADS = "1")
      if (requireNamespace("RhpcBLASctl", quietly = TRUE))
        RhpcBLASctl::blas_set_num_threads(1L)
      NULL
    })
  }

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
      cases <- list(
        no_effect = d0,
        constant_effect = make_constant(d0, d1$truth$visit_scale, delta),
        heterogeneous_effect = d1
      )

      seed <- 200001L + (B + 1L) * sim
      obj <- lapply(names(cases), function(e)
        run_condition(cases[[e]], seed, zero = (e == "no_effect"), cl = cl))
      names(obj) <- names(cases)
      saveRDS(obj, path)
    }

    all_results[[sim]] <- do.call(rbind, lapply(names(obj), function(e) {
      z <- obj[[e]]; z$effect <- e; z$sim <- sim; z
    }))

    if (sim %% 10L == 0L || sim == M)
      write.csv(do.call(rbind, all_results[seq_len(sim)]),
                file.path(OUT, "results.csv"), row.names = FALSE)
  }
}, finally = {if (!is.null(cl)) parallel::stopCluster(cl)})

results <- do.call(rbind, all_results)
saveRDS(results, file.path(OUT, "results.rds"))

# Calibration: Visit 1, Visit 2, Visit 3, Mean
visits <- c("visit1", "visit2", "visit3", "mean")
grid <- seq(0, 1, .01)
curve <- function(p) 100 * vapply(grid, function(a) mean(p <= a), 0.0)

png(file.path(OUT, "Calibration.png"), width = 1800, height = 1400, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.5, 1))

for (v in visits) {
  no <- subset(results, effect == "no_effect" & visit == v)
  co <- subset(results, effect == "constant_effect" & visit == v)

  panel <- data.frame(
    Nominal_Percent = 100 * grid,
    ATE_ZeroEffect_Percent = curve(no$ATE_p_positive),
    CATE_ConstantEffect_Percent = curve(co$CATE_p_homogeneity)
  )
  write.csv(panel, file.path(OUT, paste0("Calibration_", v, ".csv")),
            row.names = FALSE)

  plot(panel$Nominal_Percent, panel$ATE_ZeroEffect_Percent,
       type = "n", xlim = c(0, 100), ylim = c(0, 100),
       xlab = "Nominal significance level (%)",
       ylab = "Empirical rejection rate (%)",
       main = if (v == "mean") "Mean" else sub("visit", "Visit ", v))
  abline(0, 1, lty = 2, col = "gray50")
  lines(panel$Nominal_Percent, panel$ATE_ZeroEffect_Percent,
        col = "blue", lwd = 2)
  lines(panel$Nominal_Percent, panel$CATE_ConstantEffect_Percent,
        col = "red", lwd = 2)

  if (v == "visit1")
    legend("topleft", c("ATE", "CATE", "Ideal"),
           col = c("blue", "red", "gray50"), lty = c(1, 1, 2),
           lwd = c(2, 2, 1), cex = .8, bty = "n")
}
par(op)
dev.off()

# Type I error and power with exact 95% binomial confidence intervals.
rate <- function(e, v, col, label) {
  d <- subset(results, effect == e & visit == v)
  x <- d[[col]] <= LEVEL
  ci <- 100 * unname(binom.test(sum(x), length(x))$conf.int)
  r <- 100 * mean(x)
  data.frame(Test = label, Rate = r, Lower = ci[1], Upper = ci[2],
             Error_Minus = r - ci[1], Error_Plus = ci[2] - r)
}

TypeI <- do.call(rbind, lapply(visits, function(v) rbind(
  rate("no_effect", v, "ATE_p_positive", paste("ATE null -", v)),
  rate("constant_effect", v, "CATE_p_homogeneity", paste("CATE constant -", v))
)))

Power <- do.call(rbind, lapply(visits, function(v) rbind(
  rate("constant_effect", v, "ATE_p_positive", paste("ATE constant -", v)),
  rate("heterogeneous_effect", v, "ATE_p_positive",
       paste("ATE heterogeneous -", v)),
  rate("heterogeneous_effect", v, "CATE_p_homogeneity",
       paste("CATE heterogeneous -", v))
)))

write.csv(TypeI, file.path(OUT, "TypeI.csv"), row.names = FALSE)
write.csv(Power, file.path(OUT, "Power.csv"), row.names = FALSE)
cat("Complete:", M, "replications x", B, "randomizations. Output:", OUT, "\n")
