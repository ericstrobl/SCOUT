# COMPACT-like latent-confounder / noisy-proxy simulation generator.

.signed_unif <- function(n, lo=.25, hi=1)
  sample(c(-1,1),n,replace=TRUE)*runif(n,lo,hi)

make_master_data <- function(
    seed,
    n_max=2000, p_max=500,
    q=20, r=20, m=3,
    proxy_r2=.50,
    rho_item=.50, rho_time=.50,
    visit_scale=c(.50,.75,1),
    hetero_strength=.50,
    p_modifier=10) {
  
  stopifnot(length(visit_scale)==m, p_modifier<=p_max)
  set.seed(seed)
  
  # Latent clinical/confounding state.
  C <- matrix(rnorm(n_max*q),n_max,q)
  
  # Treatment mechanism.
  a <- .signed_unif(q)/sqrt(q)
  e <- plogis(drop(C%*%a))
  T <- rbinom(n_max,1,e)       # T=1 untreated; T=0 treated
  
  # Noisy observed proxies of C.
  Lambda <- matrix(rnorm(q*p_max),q,p_max)
  Lambda <- sweep(Lambda,2,sqrt(colSums(Lambda^2)),"/")
  
  X <- sqrt(proxy_r2)*(C%*%Lambda) +
    sqrt(1-proxy_r2)*matrix(rnorm(n_max*p_max),n_max,p_max)
  
  colnames(X) <- paste0("X",seq_len(p_max))
  
  # Observed effect modifier M = X b, chosen to be treatment-relevant.
  jj <- seq_len(p_modifier)
  Lm <- Lambda[,jj,drop=FALSE]
  
  Sigma_m <- proxy_r2*crossprod(Lm) +
    (1-proxy_r2)*diag(p_modifier)
  
  cov_XT <- sqrt(proxy_r2)*drop(crossprod(Lm,a))
  
  b_mod <- drop(solve(Sigma_m,cov_XT))
  b_mod <- b_mod /
    sqrt(drop(crossprod(b_mod,Sigma_m%*%b_mod)))
  
  M <- drop(X[,jj,drop=FALSE]%*%b_mod)   # Var(M)=1 in population
  
  b_full <- numeric(p_max)
  b_full[jj] <- b_mod
  
  # Baseline outcome effects of latent state.
  B <- matrix(.signed_unif(q*r),q,r)/sqrt(q)
  
  for (j in seq_len(r))
    if (sum(a*B[,j])<0)
      B[,j] <- -B[,j]
  
  # Marginal treatment effects.
  delta <- .signed_unif(r)
  while (all(delta<=0) || all(delta>=0))
    delta <- .signed_unif(r)
  
  # Heterogeneous-effect loadings.
  # Positive so every nonnegative alpha retains genuine heterogeneity.
  zeta <- runif(r,.5,1.5)
  zeta <- zeta/sqrt(mean(zeta^2))
  
  # Correlated outcome noise.
  S_item <- outer(seq_len(r),seq_len(r),
                  function(i,j) rho_item^abs(i-j))
  S_time <- outer(seq_len(m),seq_len(m),
                  function(i,j) rho_time^abs(i-j))
  
  E <- matrix(rnorm(n_max*m*r),n_max,m*r) %*%
    chol(kronecker(S_time,S_item))
  
  muC <- C%*%B
  
  Y_null <- Y_alt <- vector("list",m)
  delta_alt <- zeta_alt <- tau_alt <- vector("list",m)
  
  for (k in seq_len(m)) {
    
    cols <- ((k-1)*r+1):(k*r)
    
    delta_alt[[k]] <- visit_scale[k]*delta
    zeta_alt[[k]]  <- visit_scale[k]*zeta
    
    # True individual treatment effects:
    #
    # tau_ikj = delta_kj +
    #           hetero_strength * M_i * zeta_kj
    tau_alt[[k]] <-
      matrix(delta_alt[[k]],n_max,r,byrow=TRUE) +
      hetero_strength*(M %o% zeta_alt[[k]])
    
    Y_null[[k]] <- muC + E[,cols,drop=FALSE]
    
    Y_alt[[k]] <- Y_null[[k]] +
      sweep(tau_alt[[k]],1,T,"*")
    
    colnames(Y_null[[k]]) <-
      colnames(Y_alt[[k]]) <-
      paste0("Y",seq_len(r))
  }
  
  # E[(T-E[T|C])^2]
  eta_sd <- sqrt(sum(a^2))
  
  vt <- integrate(function(z) {
    ez <- plogis(z)
    ez*(1-ez)*dnorm(z,sd=eta_sd)
  },-Inf,Inf)$value
  
  list(
    X=X,T=T,Y_null=Y_null,Y_alt=Y_alt,
    
    truth=list(
      q=q,r=r,m=m,
      a=a,B=B,Lambda=Lambda,
      
      delta_alt=delta_alt,
      delta_null=lapply(delta_alt,\(x) rep(0,length(x))),
      
      zeta_alt=zeta_alt,
      
      tau_alt=tau_alt,
      
      modifier=M,
      modifier_beta=b_full,
      p_modifier=p_modifier,
      hetero_strength=hetero_strength,
      
      Sigma_item=S_item,
      treatment_resid_var=vt,
      proxy_r2=proxy_r2,
      visit_scale=visit_scale
    )
  )
}


get_sim_data <- function(master,n,p,
                         effect=c("alternative","null")) {
  
  effect <- match.arg(effect)
  
  if (n>nrow(master$X) || p>ncol(master$X))
    stop("Requested n or p exceeds the master dataset.")
  
  if (p<master$truth$p_modifier)
    stop("p must be >= p_modifier.")
  
  Y <- if (effect=="alternative")
    master$Y_alt
  else
    master$Y_null
  
  truth <- master$truth
  
  # Keep truth corresponding to observed subjects/covariates.
  truth$Lambda <- truth$Lambda[,seq_len(p),drop=FALSE]
  truth$modifier_beta <- truth$modifier_beta[seq_len(p)]
  truth$modifier <- truth$modifier[seq_len(n)]
  
  if (effect=="alternative") {
    
    truth$delta <- truth$delta_alt
    truth$zeta <- truth$zeta_alt
    
    truth$tau <- lapply(
      truth$tau_alt,
      \(x) x[seq_len(n),,drop=FALSE]
    )
    
  } else {
    
    truth$delta <- truth$delta_null
    truth$zeta <- lapply(
      truth$zeta_alt,
      \(x) rep(0,length(x))
    )
    
    truth$tau <- lapply(
      seq_len(truth$m),
      \(k) matrix(0,n,truth$r)
    )
  }
  
  truth$effect <- effect
  
  list(
    X=master$X[seq_len(n),seq_len(p),drop=FALSE],
    T=master$T[seq_len(n)],
    Y_list=lapply(
      Y,
      \(x) x[seq_len(n),,drop=FALSE]
    ),
    truth=truth
  )
}