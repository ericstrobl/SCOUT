make_dgp <- function(seed,
                     q=20, r=20, p_max=500,
                     proxy_r2=.50,
                     visit_scale=c(.50,.75,1),
                     p_modifier=10) {
  
  set.seed(seed)
  
  # The treatment index is a latent direction, not a function of observed X.
  a <- .signed_unif(q) / sqrt(q)
  
  Lambda <- matrix(rnorm(q*p_max), q, p_max)
  Lambda <- sweep(Lambda, 2, sqrt(colSums(Lambda^2)), "/")
  
  # Observed effect modifier M = X[, 1:p_modifier] %*% b_mod, chosen to be
  # correlated with the latent treatment index a'C (not a balancing score).
  jj <- seq_len(p_modifier)
  Lm <- Lambda[, jj, drop=FALSE]
  Sigma_m <- proxy_r2 * crossprod(Lm) +
    (1 - proxy_r2) * diag(p_modifier)
  cov_XU <- sqrt(proxy_r2) * drop(crossprod(Lm, a))
  b_mod <- drop(solve(Sigma_m, cov_XU))
  b_mod <- b_mod / sqrt(drop(crossprod(b_mod, Sigma_m %*% b_mod)))
  
  B <- matrix(.signed_unif(q*r), q, r) / sqrt(q)
  for (j in seq_len(r))
    if (sum(a*B[,j]) < 0) B[,j] <- -B[,j]
  
  delta <- .signed_unif(r)
  while (all(delta<=0) || all(delta>=0))
    delta <- .signed_unif(r)
  
  zeta <- runif(r,.5,1.5)
  zeta <- zeta / sqrt(mean(zeta^2))
  
  list(
    q=q, r=r, p_max=p_max,
    a=a, Lambda=Lambda, B=B,
    delta=delta, zeta=zeta,
    b_mod=b_mod,
    proxy_r2=proxy_r2,
    visit_scale=visit_scale,
    p_modifier=p_modifier
  )
}


sample_from_dgp <- function(
    dgp, seed, n_max=2000, m=3,
    rho_item=.50, rho_time=.50,
    hetero_strength=.50) {
  
  set.seed(seed)
  
  q <- dgp$q
  r <- dgp$r
  p_max <- dgp$p_max
  
  C <- matrix(rnorm(n_max*q), n_max, q)
  
  # Noisy observed proxies of the latent state (same model as main simulations).
  X <- sqrt(dgp$proxy_r2) * (C %*% dgp$Lambda) +
    sqrt(1-dgp$proxy_r2) *
    matrix(rnorm(n_max*p_max), n_max, p_max)
  
  colnames(X) <- paste0("X",seq_len(p_max))
  
  # Treatment is logistic in the unobserved a'C, not in X or any observed
  # score. Residual uncertainty about C remains after conditioning on X.
  e <- plogis(drop(C %*% dgp$a))
  T <- rbinom(n_max, 1, e)
  
  # Observed proxy-based effect modifier, correlated with latent treatment risk.
  jj <- seq_len(dgp$p_modifier)
  M <- drop(X[, jj, drop=FALSE] %*% dgp$b_mod)
  
  S_item <- outer(
    seq_len(r),seq_len(r),
    function(i,j) rho_item^abs(i-j)
  )
  
  S_time <- outer(
    seq_len(m),seq_len(m),
    function(i,j) rho_time^abs(i-j)
  )
  
  E <- matrix(rnorm(n_max*m*r),n_max,m*r) %*%
    chol(kronecker(S_time,S_item))
  
  muC <- C %*% dgp$B
  
  Y_null <- Y_alt <- vector("list",m)
  delta_alt <- zeta_alt <- tau_alt <- vector("list",m)
  
  for (k in seq_len(m)) {
    
    cols <- ((k-1)*r+1):(k*r)
    
    delta_alt[[k]] <- dgp$visit_scale[k] * dgp$delta
    zeta_alt[[k]]  <- dgp$visit_scale[k] * dgp$zeta
    
    tau_alt[[k]] <-
      matrix(delta_alt[[k]],n_max,r,byrow=TRUE) +
      hetero_strength * (M %o% zeta_alt[[k]])
    
    Y_null[[k]] <- muC + E[,cols,drop=FALSE]
    Y_alt[[k]] <- Y_null[[k]] +
      sweep(tau_alt[[k]],1,T,"*")
  }
  
  list(
    X=X, T=T, Y_null=Y_null, Y_alt=Y_alt,
    truth=list(
      q=q, r=r, m=m,
      a=dgp$a,
      B=dgp$B,
      Lambda=dgp$Lambda,
      delta_alt=delta_alt,
      delta_null=lapply(delta_alt,\(x) rep(0,length(x))),
      zeta_alt=zeta_alt,
      tau_alt=tau_alt,
      modifier=M,
      p_modifier=dgp$p_modifier,
      hetero_strength=hetero_strength,
      proxy_r2=dgp$proxy_r2,
      visit_scale=dgp$visit_scale
    )
  )
}

