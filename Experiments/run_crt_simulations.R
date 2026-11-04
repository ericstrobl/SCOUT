
# run_crt_simulations.R
# Monte Carlo evaluation of SCOUT's ATE and CATE CRTs.
# Uses the functions in R/ and saves results after each replication.

SRC_DIR <- "R"
for (file in c("sim_data.R", "randomization_dgp.R",
               "SCOUT.R", "SCOUT_CRT.R"))
  source(file.path(SRC_DIR, file))

# Simulation settings
M <- 500L
B <- 500L
n <- 1000L
p <- 200L
DGP_SEED <- 12345L
VERSION <- 25L
CRT_CORES <- 4L
CRT_RESELECT_ETA <- TRUE
LEVEL <- 0.05

OUT <- sprintf("scout_CRT_M%d_B%d_n%d_p%d_dgp%d_v%d",
               M, B, n, p, DGP_SEED, VERSION)
dir.create(OUT, showWarnings = FALSE)

config <- list(M=M, B=B, n=n, p=p, DGP_SEED=DGP_SEED,
               VERSION=VERSION, CRT_RESELECT_ETA=CRT_RESELECT_ETA)
config_file <- file.path(OUT, "config.rds")

if (file.exists(config_file)) {
  if (!identical(readRDS(config_file), config))
    stop("Cached settings differ; change VERSION.")
} else {
  saveRDS(config, config_file)
}

# Constant outcome-specific effects sampled independently from U(-1,1).
# Effects may vary across outcomes and visits, but not individuals.
make_constant <- function(dat, seed) {
  set.seed(seed)
  out <- dat
  nn <- nrow(dat$X)
  r <- ncol(dat$Y_list[[1]])
  delta <- runif(r, -1, 1)
  
  for (k in seq_along(out$Y_list)) {
    dk <- dat$truth$visit_scale[k] * delta
    out$Y_list[[k]] <- dat$Y_list[[k]] + tcrossprod(dat$T, dk)
    out$truth$tau[[k]] <- matrix(rep(dk, each=nn), nn, r)
  }
  
  out$truth$effect <- "constant"
  out
}

# Fit SCOUT, then obtain unadjusted CRT p-values.
run_condition <- function(dat, seed) {
  fit <- SCOUT(dat$X, dat$T, dat$Y_list, seed=1L)
  
  crt <- SCOUT_CRT(
    fit, dat$X, dat$T, dat$Y_list,
    B=B, seed=seed, cores=CRT_CORES,
    reselect_eta=CRT_RESELECT_ETA,
    fit_seed=1L
  )
  
  crt$results
}

# Fixed DGP parameters; independent samples across replications.
dgp <- make_dgp(DGP_SEED)
saveRDS(dgp, file.path(OUT, "fixed_dgp.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT, "sessionInfo.txt"))

all_results <- vector("list", M)

for (sim in seq_len(M)) {
  message("Replication ", sim, "/", M)
  path <- file.path(OUT, sprintf("rep_%03d.rds", sim))
  
  if (file.exists(path)) {
    z <- readRDS(path)
  } else {
    master <- sample_from_dgp(dgp, seed=10001L+sim, n_max=n)
    d0 <- get_sim_data(master, n, p, "null")
    d1 <- get_sim_data(master, n, p, "alternative")
    
    conditions <- list(
      no_effect = d0,
      constant_effect = make_constant(d0, 300001L+sim),
      heterogeneous_effect = d1
    )
    
    z <- do.call(rbind, lapply(names(conditions), function(effect) {
      result <- run_condition(
        conditions[[effect]],
        200001L + (B + 1L)*sim
      )
      result$effect <- effect
      result$sim <- sim
      result
    }))
    
    rownames(z) <- NULL
    saveRDS(z, path)
  }
  
  all_results[[sim]] <- z
  
  if (sim %% 10L == 0L || sim == M)
    write.csv(do.call(rbind, all_results[seq_len(sim)]),
              file.path(OUT, "results.csv"), row.names=FALSE)
}

results <- do.call(rbind, all_results)
rownames(results) <- NULL
saveRDS(results, file.path(OUT, "results.rds"))
write.csv(results, file.path(OUT, "results.csv"), row.names=FALSE)

# Rejection rates and exact binomial 95% confidence intervals.
summarize_group <- function(d) {
  do.call(rbind, lapply(c("ATE", "CATE"), function(test) {
    reject <- d[[paste0(test, "_p")]] <= LEVEL
    ci <- binom.test(sum(reject), length(reject))$conf.int
    
    data.frame(
      effect=d$effect[1], visit=d$visit[1], test=test,
      M_completed=nrow(d),
      Rate=100*mean(reject),
      Lower=100*ci[1], Upper=100*ci[2]
    )
  }))
}

summary_results <- do.call(rbind, lapply(
  split(results, interaction(results$effect, results$visit)),
  summarize_group
))
rownames(summary_results) <- NULL

write.csv(summary_results, file.path(OUT, "summary.csv"),
          row.names=FALSE)

# Calibration curves for the two null hypotheses.
grid <- seq(.01, 1, .01)

get_p <- function(effect, column)
  results[results$effect == effect &
            results$visit == "mean_visit_Q", column]

curve <- function(pvals)
  100*vapply(grid, function(a) mean(pvals <= a), numeric(1))

panelA <- data.frame(
  Nominal_Percent=100*grid,
  ATE_ZeroEffect_Percent=curve(get_p("no_effect", "ATE_p")),
  CATE_ZeroEffect_Percent=curve(get_p("no_effect", "CATE_p")),
  CATE_ConstantEffect_Percent=curve(get_p("constant_effect", "CATE_p"))
)

write.csv(panelA, file.path(OUT, "Prism_Panel_A_Calibration.csv"),
          row.names=FALSE)

# Type I error: ATE under zero effects; CATE under homogeneous effects.
panelB <- subset(summary_results,
                 (effect == "no_effect" & test == "ATE") |
                   (effect %in% c("no_effect", "constant_effect") & test == "CATE")
)

# Power: ATE under nonzero effects; CATE under heterogeneous effects.
panelC <- subset(summary_results,
                 (effect %in% c("constant_effect", "heterogeneous_effect") &
                    test == "ATE") |
                   (effect == "heterogeneous_effect" & test == "CATE")
)

write.csv(panelB, file.path(OUT, "Prism_Panel_B_TypeI.csv"),
          row.names=FALSE)
write.csv(panelC, file.path(OUT, "Prism_Panel_C_Power.csv"),
          row.names=FALSE)

print(summary_results, row.names=FALSE)
cat("\nComplete:", M, "replications,", B,
    "randomizations per condition.\nOutput:", OUT, "\n")
