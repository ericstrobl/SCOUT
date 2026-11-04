
# SCOUT

SCOUT (Score-based Causal Outcome Learning Using Treatment Signals) is an algorithm for learning treatment-responsive composite outcomes while adjusting for confounding in observational biomedical data. SCOUT addresses settings with numerous candidate outcomes and high-dimensional, noisy baseline variables that serve as proxies for potentially unobserved confounders. The algorithm first constructs a common adjustment representation for learning nonnegative outcome weights across multiple longitudinal visits. After learning the composite outcome, SCOUT reduces the adjustment representation to at most two dimensions using COMPACT, enabling estimation of average treatment effects (ATEs) and conditional average treatment effects (CATEs) within the learned score space. SCOUT also provides conditional randomization tests (CRTs) for detecting average treatment effects and treatment-effect heterogeneity while accounting for adaptive outcome learning.

Code was tested in R version 4.3.1.

# Installation

> install.packages("remotes")

> remotes::install_github("ericstrobl/SCOUT")

> library(SCOUT)

# Sample Run

Generate synthetic observational data with latent confounding, noisy baseline proxies, and heterogeneous treatment effects:

```r
n <- 300
p <- 30
r <- 5
m <- 3

C <- matrix(rnorm(n * 3), n, 3)
X <- cbind(C + matrix(rnorm(n * 3), n, 3),
           matrix(rnorm(n * (p - 3)), n, p - 3))
T <- rbinom(n, 1, plogis(0.6*C[,1] - 0.4*C[,2]))

Y_list <- lapply(seq_len(m), function(k) {
  Y <- sapply(seq_len(r), function(j) {
    0.4*C[,1] + 0.2*C[,2] +
      T*(0.3 + 0.1*j + 1*C[,1])*k/m + rnorm(n)
  })
  colnames(Y) <- paste0("Y", seq_len(r))
  Y
})
```

Run SCOUT with cross-validation:

> fit <- SCOUT(X, T, Y_list)

View the learned outcome weights:

> fit$alpha

View the two-dimensional adjustment scores at the first visit:

> fit$final[[1]]$S_visit

Estimated ATEs across visits:

> sapply(fit$final, function(x) x$ATE)

Estimated CATEs at the first visit:

> fit$final[[1]]$CATE

# Conditional Randomization Tests

Run CRT inference with 500 randomizations:

> crt <- SCOUT_CRT(fit, X, T, Y_list, B = 500, cores = 4)

View the results:

> crt$results

The results contain the following columns:

- `visit`: Follow-up visit or pooled analysis (`mean_visit_Q`).
- `n_visit`: Number of observed patients at each visit.
- `ATE`: Estimated average treatment effect.
- `ATE_p`: One-sided p-value for a positive average treatment effect.
- `CATE_Q`: Fourth-order treatment-effect heterogeneity statistic.
- `CATE_p`: P-value for treatment-effect heterogeneity under the constant-effect null.

The `mean_visit_Q` row tests the mean ATE and the mean of the visit-specific fourth-order heterogeneity statistics across visits.

# Incomplete Longitudinal Data

SCOUT also supports different patients observed at different visits. The `visit_ids` argument specifies which baseline individuals correspond to each outcome matrix.

For example:

```r
set.seed(123)

visit_ids <- list(
  seq_len(n),
  sort(sample(seq_len(n), 250)),
  sort(sample(seq_len(n), 200))
)

Y_incomplete <- Map(
  function(Y, ids) Y[ids, , drop = FALSE],
  Y_list, visit_ids
)
```

Fit SCOUT and run CRT inference:

> fit <- SCOUT(X, T, Y_incomplete, visit_ids = visit_ids)

> crt <- SCOUT_CRT(fit, X, T, Y_incomplete, B = 500, cores = 4)

> crt$results

The CRT retains a single treatment assignment per patient across visits while allowing the observed sample to differ between follow-ups.

# Reproducing the Simulation Results

The `Experiments` folder contains the scripts used to evaluate SCOUT against alternative methods and assess CRT calibration and statistical power.

- `run_simulations.R`: Synthetic experiments comparing ATE estimation, CATE estimation, and computational runtime across sample sizes and baseline dimensions.
- `run_crt_simulations.R`: Monte Carlo experiments evaluating CRT Type I error, calibration, and statistical power under zero, constant, and heterogeneous treatment effects.
