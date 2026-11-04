# Exact metrics for the main point-performance simulation.

.simplex_proj <- function(x) {
  u <- sort(x, decreasing = TRUE)
  z <- cumsum(u) - 1
  j <- max(which(u - z / seq_along(u) > 0))
  pmax(x - z[j] / j, 0)
}

true_J0 <- function(alpha, truth) {
  if (any(!is.finite(alpha))) return(NA_real_)
  V <- truth$treatment_resid_var
  S <- truth$Sigma_item
  
  sum(vapply(truth$delta, function(dlt) {
    d <- sum(alpha * dlt)
    v <- drop(crossprod(alpha, S %*% alpha))
    den <- V * d^2 + v
    if (den <= 0) 0 else sqrt(V) * d / sqrt(den)
  }, 0.0))
}

.J0_grad <- function(alpha, truth) {
  V <- truth$treatment_resid_var
  S <- truth$Sigma_item
  Sa <- drop(S %*% alpha)
  v <- sum(alpha * Sa)
  g <- numeric(length(alpha))
  
  for (dlt in truth$delta) {
    d <- sum(alpha * dlt)
    D <- V * d^2 + v
    if (D > 1e-14)
      g <- g + sqrt(V) * (v * dlt - d * Sa) / D^(3/2)
  }
  g
}

oracle_alpha <- function(truth, starts = 10, seed = 1,
                         maxit = 1000, tol = 1e-9) {
  r <- truth$r
  
  set.seed(seed)
  A <- list(rep(1 / r, r))
  
  vertex <- vapply(seq_len(r), function(j) {
    a <- numeric(r); a[j] <- 1
    true_J0(a, truth)
  }, 0.0)
  
  a <- numeric(r); a[which.max(vertex)] <- 1
  A[[2]] <- a
  
  if (starts > 2)
    for (j in 3:starts) {
      a <- rexp(r)
      A[[j]] <- a / sum(a)
    }
  
  one <- function(a) {
    a <- .simplex_proj(a)
    q <- true_J0(a, truth)
    
    for (it in seq_len(maxit)) {
      g <- .J0_grad(a, truth)
      step <- 1
      
      for (ls in 1:50) {
        b <- .simplex_proj(a + step * g)
        qb <- true_J0(b, truth)
        if (is.finite(qb) && qb >= q - 1e-12) break
        step <- step / 2
      }
      
      if (sum(abs(b - a)) < tol || abs(qb - q) < tol) {
        a <- b
        q <- qb
        break
      }
      
      a <- b
      q <- qb
    }
    
    list(alpha = a, objective = q)
  }
  
  z <- lapply(A, one)
  z[[which.max(vapply(z, `[[`, 0.0, "objective"))]]$alpha
}

metric_reference <- function(truth, seed = 1) {
  af <- oracle_alpha(truth, seed = seed)
  
  # Population covariance of the questionnaire within a fixed treatment arm:
  # Var(C B + epsilon). Additive treatment shifts do not affect covariance.
  Sph <- crossprod(truth$B) + truth$Sigma_item
  
  list(alpha_full = af, Sigma_pheno = Sph)
}

phenotype_similarity <- function(alpha, ref) {
  b <- ref$alpha_full
  
  if (any(!is.finite(alpha)) || any(!is.finite(b)))
    return(NA_real_)
  
  S <- ref$Sigma_pheno
  num <- drop(crossprod(alpha, S %*% b))
  den <- sqrt(
    drop(crossprod(alpha, S %*% alpha)) *
      drop(crossprod(b, S %*% b))
  )
  
  if (den <= 0) NA_real_ else num / den
}

true_effect <- function(alpha, truth, visit = truth$m) {
  if (any(!is.finite(alpha))) return(NA_real_)
  sum(alpha * truth$delta[[visit]])
}

metric_row <- function(sim, n, p, algorithm,
                       alpha, estimate, runtime_sec,
                       truth, ref) {
  j0 <- true_J0(alpha, truth)
  target <- true_effect(alpha, truth)
  
  data.frame(
    sim = sim,
    n = n,
    p = p,
    algorithm = algorithm,
    J0 = j0,
    signal_positive = is.finite(j0) && j0 > 0,
    phenotype_similarity = phenotype_similarity(alpha, ref),
    estimate = estimate,
    true_effect = target,
    error = if (is.finite(estimate) && is.finite(target))
      estimate - target else NA_real_,
    squared_error = if (is.finite(estimate) && is.finite(target))
      (estimate - target)^2 else NA_real_,
    runtime_sec = runtime_sec,
    status = "ok",
    message = "",
    stringsAsFactors = FALSE
  )
}

metric_failure <- function(sim, n, p, algorithm, msg) {
  data.frame(
    sim = sim,
    n = n,
    p = p,
    algorithm = algorithm,
    J0 = NA_real_,
    signal_positive = NA,
    phenotype_similarity = NA_real_,
    estimate = NA_real_,
    true_effect = NA_real_,
    error = NA_real_,
    squared_error = NA_real_,
    runtime_sec = NA_real_,
    status = "error",
    message = as.character(msg),
    stringsAsFactors = FALSE
  )
}

summarize_metrics <- function(results) {
  key <- interaction(
    results$algorithm, results$n, results$p,
    drop = TRUE
  )
  groups <- split(results, key)
  
  avg <- function(x)
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  
  do.call(rbind, lapply(groups, function(d) {
    data.frame(
      algorithm = d$algorithm[1],
      n = d$n[1],
      p = d$p[1],
      causal_signal = avg(d$J0),
      signal_preservation = avg(d$signal_positive),
      phenotype_similarity = avg(d$phenotype_similarity),
      bias = avg(d$error),
      RMSE = sqrt(avg(d$squared_error)),
      runtime_sec = avg(d$runtime_sec),
      failures = sum(d$status != "ok")
    )
  }))
}

# Cache point-fit quantities and exact truth so additional metrics can be
# derived later without rerunning any algorithm.
make_metric_cache <- function(
    sim, n, p, algorithm, data_seed, algorithm_seed,
    truth, ref, point, diagnostics = NULL) {
  
  list(
    meta = list(
      sim = sim,
      n = n,
      p = p,
      algorithm = algorithm,
      data_seed = data_seed,
      algorithm_seed = algorithm_seed
    ),
    truth = truth,
    reference = ref,
    point = point,
    diagnostics = diagnostics
  )
}
