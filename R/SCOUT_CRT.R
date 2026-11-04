# SCOUT_CRT.R — CRT inference for an existing SCOUT() fit.

.scoutcrt_ids <- function(X, T, Y, ids) {
  n <- nrow(X)
  if (is.null(ids)) {
    if (any(vapply(Y, nrow, 0L) != n))
      stop("Supply visit_ids for incomplete visits.")
    ids <- rep(list(seq_len(n)), length(Y))
  }
  if (!is.list(ids) || length(ids) != length(Y))
    stop("visit_ids must have one index vector per visit.")
  for (k in seq_along(Y)) {
    ii <- ids[[k]]
    if (!is.numeric(ii) || anyNA(ii) || any(!is.finite(ii)) ||
        any(ii != as.integer(ii)) || any(ii < 1 | ii > n) ||
        anyDuplicated(ii) || length(ii) != nrow(Y[[k]]))
      stop("Invalid visit_ids at visit ", k)
    if (length(ii) < 3L || length(unique(T[ii])) != 2L)
      stop("Visit ", k, " requires at least 3 subjects and both treatments.")
    ids[[k]] <- as.integer(ii)
  }
  ids
}

# Per-visit treatment data, independent of outcomes; reused for both nulls.
.scoutcrt_prep_xt <- function(X, T, ids) {
  x <- X[ids, , drop = FALSE]
  t <- T[ids]
  n <- length(ids)
  xm <- colMeans(x)
  x <- sweep(x, 2, xm, "-")
  tc <- t - mean(t)
  tt <- sum(tc^2)
  if (tt <= 0) stop("Treatment does not vary within visit.")
  xt <- crossprod(x, tc)
  G <- crossprod(x)/n - tcrossprod(drop(xt))/(2*n*tt)
  E <- eigen((G + t(G))/2, symmetric = TRUE)
  list(X = x, T = t, ids = ids, n = n, p = ncol(X), xm = xm,
       tc = tc, tt = tt, U = E$vectors,
       raw_ev = E$values, ev = pmax(E$values, 0),
       uT = crossprod(E$vectors, xt/n))
}

.scoutcrt_prep_y <- function(P, Y) {
  yc <- sweep(Y, 2, colMeans(Y), "-")
  ty <- crossprod(P$tc, yc)
  zy <- crossprod(P$U, crossprod(P$X, yc)/P$n)
  list(Y = Y, TY = ty, ZY = zy,
       V = zy - P$uT %*% (ty/P$tt), YY = crossprod(yc))
}

# Whiten scores using small matrices in the eigenspace of G.
.scoutcrt_basis <- function(P, Q, eta, with_B = FALSE) {
  D <- sweep(Q, 1, P$ev + eta, "/")
  C <- crossprod(D, sweep(D, 1, P$raw_ev, "*"))
  u <- drop(crossprod(D, P$uT))
  C <- C + (P$n/(2*P$tt))*tcrossprod(u)
  E <- eigen((C + t(C))/2, symmetric = TRUE)
  mx <- max(E$values)
  keep <- if (is.finite(mx) && mx > 0) E$values > 1e-10*mx else
    rep(FALSE, length(E$values))
  if (any(keep)) {
    W <- sweep(E$vectors[, keep, drop = FALSE],
               2, sqrt(E$values[keep]), "/")
    R <- D %*% W
    gram <- crossprod(W, C %*% W)
  } else {
    R <- matrix(0, P$p, 0)
    gram <- matrix(0, 0, 0)
  }
  list(R = R, gram = gram,
       B = if (with_B) P$U %*% R else NULL)
}

.scoutcrt_nnls <- function(H, cvec, tt) {
  f <- function(a) tt - 2*sum(a*cvec) + sum(a*drop(H %*% a))
  g <- function(a) 2*(drop(H %*% a) - cvec)
  out <- stats::optim(rep(0, length(cvec)), f, g, method = "L-BFGS-B",
                      lower = rep(0, length(cvec)),
                      control = list(maxit = 500, pgtol = 1e-8))
  a <- pmax(out$par, 0)
  if (any(!is.finite(a))) stop("NNLS produced nonfinite weights.")
  a
}

# Accelerated SCOUT refit: cross-products instead of residualized n-by-r Y.
.scoutcrt_fit <- function(P, O, eta, effects = TRUE) {
  K <- length(P)
  r <- ncol(O[[1]]$Y)
  H <- matrix(0, r, r)
  cvec <- numeric(r)
  tt <- 0
  for (k in seq_len(K)) {
    pk <- P[[k]]
    ok <- O[[k]]
    z <- .scoutcrt_basis(pk, cbind(pk$uT, ok$V), eta)
    if (!ncol(z$R)) {
      Hk <- ok$YY
      ck <- drop(ok$TY)
      tk <- pk$tt
    } else {
      sy <- crossprod(z$R, ok$ZY)
      st <- drop(crossprod(z$R, pk$uT))
      G <- z$gram
      Hk <- ok$YY - 2*pk$n*crossprod(sy) +
        pk$n*crossprod(sy, G %*% sy)
      ck <- drop(ok$TY) - 2*pk$n*drop(crossprod(sy, st)) +
        pk$n*drop(crossprod(sy, G %*% st))
      tk <- pk$tt - 2*pk$n*sum(st^2) +
        pk$n*drop(crossprod(st, G %*% st))
    }
    # Equal weight per visit, regardless of n_k.
    H <- H + Hk/pk$n
    cvec <- cvec + ck/pk$n
    tt <- tt + tk/pk$n
  }
  H <- (H + t(H))/2
  if (any(!is.finite(H)) || any(!is.finite(cvec)) ||
      !is.finite(tt) || tt <= .Machine$double.eps)
    stop("Invalid residualized NNLS cross-products.")
  a <- .scoutcrt_nnls(H, cvec, tt)
  if (sum(a) > 0) a <- a/sum(a)
  
  final <- lapply(seq_len(K), function(k) {
    pk <- P[[k]]
    ok <- O[[k]]
    B <- .scoutcrt_basis(pk, cbind(pk$uT, ok$V %*% a),
                         eta, with_B = TRUE)$B
    if (!effects) return(list(B = B))
    # Only visit-observed scores are needed for CRT statistics.
    S <- pk$X %*% B
    y <- drop(ok$Y %*% a)
    q <- ncol(S)
    coef <- stats::lm.fit(cbind(1, pk$T, S, S*pk$T), y)$coefficients
    if (any(!is.finite(coef))) stop("Rank-deficient final OLS.")
    cate <- if (q) drop(coef[2] + S %*% tail(coef, q)) else
      rep(coef[2], pk$n)
    list(B = B, ATE = mean(cate), CATE = cate)
  })
  list(alpha = a, final = final, eta = eta)
}

.scoutcrt_resid <- function(A, Z) stats::lm.fit(cbind(1, Z), A)$residuals
.scoutcrt_R2 <- function(y, Z) {
  Z <- as.matrix(Z)
  d <- sum((y - mean(y))^2)
  if (!ncol(Z) || d <= .Machine$double.eps) return(0)
  max(0, 1 - sum(.scoutcrt_resid(y, Z)^2)/d)
}

# Joint refits share the treatment eigensystems at each visit and in CV folds.
.scoutcrt_refit_pair <- function(X, T, Ysets, ids, eta = NULL,
                                 eta_grid = 10^seq(-5, 1, length.out = 12),
                                 nfolds = 5L, seed = 1L) {
  nd <- length(Ysets)
  K <- length(ids)
  if (!nd || any(vapply(Ysets, length, 0L) != K))
    stop("Invalid outcome sets.")
  
  mkP <- function(index) {
    out <- vector("list", length(index))
    for (k in seq_along(index)) {
      prev <- if (k > 1L) which(vapply(index[seq_len(k-1L)],
                                       identical, logical(1), index[[k]]))
      else integer(0)
      out[[k]] <- if (length(prev)) out[[prev[1L]]] else
        .scoutcrt_prep_xt(X, T, index[[k]])
    }
    out
  }
  mkO <- function(P, Y) Map(.scoutcrt_prep_y, P, Y)
  
  if (is.null(eta)) {
    n <- length(T)
    if (any(table(factor(T, levels = 0:1)) < nfolds) ||
        any(vapply(ids, function(ii)
          any(table(factor(T[ii], levels = 0:1)) < nfolds), logical(1))))
      stop("Too few subjects per treatment group for CV.")
    set.seed(seed)
    valid <- FALSE
    for (attempt in seq_len(200L)) {
      fold <- integer(n)
      for (g in 0:1) {
        ii <- which(T == g)
        fold[ii] <- sample(rep(seq_len(nfolds), length.out = length(ii)))
      }
      valid <- all(vapply(ids, function(ii) {
        all(vapply(seq_len(nfolds), function(f) {
          tr <- ii[fold[ii] != f]
          va <- ii[fold[ii] == f]
          length(va) >= 2L && length(unique(T[tr])) == 2L &&
            length(unique(T[va])) == 2L
        }, logical(1)))
      }, logical(1)))
      if (valid) break
    }
    if (!valid) stop("Cannot construct valid subject-level CV folds.")
    
    CV <- lapply(seq_len(nd), function(s)
      matrix(NA_real_, nfolds, length(eta_grid)))
    for (f in seq_len(nfolds)) {
      train_rows <- lapply(ids, function(ii) which(fold[ii] != f))
      test_rows <- lapply(ids, function(ii) which(fold[ii] == f))
      train_ids <- Map(function(ii, rr) ii[rr], ids, train_rows)
      test_ids <- Map(function(ii, rr) ii[rr], ids, test_rows)
      P <- mkP(train_ids)
      trainY <- lapply(Ysets, function(Y)
        Map(function(y, rr) y[rr, , drop = FALSE], Y, train_rows))
      testY <- lapply(Ysets, function(Y)
        Map(function(y, rr) y[rr, , drop = FALSE], Y, test_rows))
      Osets <- lapply(trainY, function(Y) mkO(P, Y))
      
      for (j in seq_along(eta_grid)) for (s in seq_len(nd)) {
        f0 <- .scoutcrt_fit(P, Osets[[s]], eta_grid[j], effects = FALSE)
        CV[[s]][f, j] <- mean(vapply(seq_len(K), function(k) {
          S <- sweep(X[test_ids[[k]], , drop = FALSE],
                     2, P[[k]]$xm, "-") %*% f0$final[[k]]$B
          if (!ncol(S)) return(0)
          y <- drop(testY[[s]][[k]] %*% f0$alpha)
          t <- T[test_ids[[k]]]
          sqrt(.scoutcrt_R2(t, S) *
                 .scoutcrt_R2(.scoutcrt_resid(y, t),
                              .scoutcrt_resid(S, t)))
        }, 0.0))
      }
    }
    etas <- vapply(CV, function(cv)
      eta_grid[which.max(colMeans(cv))], 0.0)
  } else {
    etas <- rep(eta, nd)
  }
  
  P <- mkP(ids)
  lapply(seq_len(nd), function(s)
    .scoutcrt_fit(P, mkO(P, Ysets[[s]]), etas[s]))
}

# Constant item-specific treatment shifts for the homogeneity null.
.scoutcrt_item_effects <- function(X, T, Y, ids, eta) {
  K <- length(ids)
  r <- ncol(Y[[1]])
  delta <- matrix(NA_real_, K, r)
  for (k in seq_len(K)) {
    P <- .scoutcrt_prep_xt(X, T, ids[[k]])
    O <- .scoutcrt_prep_y(P, Y[[k]])
    D <- cbind(P$uT, O$V)
    W <- P$U %*% sweep(D, 1, P$ev + eta, "/")
    S <- P$X %*% W
    for (j in seq_len(r))
      delta[k, j] <- stats::lm.fit(
        cbind(1, P$T, S[, 1], S[, j+1]), Y[[k]][, j]
      )$coefficients[2]
  }
  if (any(!is.finite(delta)))
    stop("Could not estimate constant item-specific treatment effects.")
  delta
}

.scoutcrt_propensity <- function(T, S, clip, supplied) {
  if (is.null(supplied)) {
    Z <- if (ncol(S)) cbind(1, S) else matrix(1, length(T), 1)
    e <- suppressWarnings(stats::glm.fit(
      Z, T, family = stats::binomial()))$fitted.values
    source <- "estimated"
  } else {
    e <- supplied
    source <- "supplied"
    if (length(e) != length(T) || any(!is.finite(e)) ||
        any(e < 0 | e > 1)) stop("Invalid supplied propensity.")
  }
  if (length(e) != length(T) || any(!is.finite(e)))
    stop("Nonfinite estimated propensity scores.")
  e <- pmin(pmax(as.numeric(e), clip), 1-clip)
  list(p = e, source = source, min = min(e), max = max(e))
}

.scoutcrt_stats <- function(f) {
  ate <- vapply(f$final, `[[`, 0.0, "ATE")
  q <- vapply(f$final, function(v) {
    d <- v$CATE - mean(v$CATE)
    sqrt(mean(d^4))
  }, 0.0)
  if (any(!is.finite(c(ate, q)))) stop("Nonfinite CRT statistics.")
  list(ATE = c(ate, mean(ate)), Q = c(q, mean(q)))
}

.scoutcrt_draw <- function(b, s) {
  set.seed(as.integer((as.numeric(s$seed) + 100000*(s$K+1L) + b) %%
                        .Machine$integer.max))
  ts <- stats::rbinom(s$n, 1L, s$e)
  out <- list(ATE = rep(NA_real_, s$J), Q = rep(NA_real_, s$J))
  if (any(vapply(s$ids, function(ii)
    min(tabulate(ts[ii]+1L, nbins = 2L)) < 5L, logical(1))))
    return(out)
  
  # Treatment-change transformation under estimated constant item effects.
  Yh <- lapply(seq_len(s$K), function(k) {
    ii <- s$ids[[k]]
    s$Y[[k]] + tcrossprod(ts[ii] - s$T[ii], s$delta[k, ])
  })
  fit_pair <- function(ys) .scoutcrt_refit_pair(
    s$X, ts, ys, s$ids, eta = s$eta,
    eta_grid = s$eta_grid, nfolds = s$nfolds, seed = s$fit_seed)
  fits <- tryCatch(fit_pair(list(s$Y, Yh)), error = function(e) NULL)
  if (is.null(fits)) {
    one <- function(y) tryCatch(fit_pair(list(y))[[1]],
                                error = function(e) NULL)
    fits <- list(one(s$Y), one(Yh))
  }
  if (!is.null(fits[[1]]))
    out$ATE <- tryCatch(.scoutcrt_stats(fits[[1]])$ATE,
                        error = function(e) rep(NA_real_, s$J))
  if (!is.null(fits[[2]]))
    out$Q <- tryCatch(.scoutcrt_stats(fits[[2]])$Q,
                      error = function(e) rep(NA_real_, s$J))
  out
}

.scoutcrt_worker <- function(b) .scoutcrt_draw(b, .scoutcrt_state)

# Public entry point. All p-values, including pooled p-values, are unadjusted.
SCOUT_CRT <- function(fit, X, T, Y_list, visit_ids = fit$visit_ids,
                      B = 500L, seed = 2026L, cores = 1L,
                      reselect_eta = !is.null(fit$eta_cv),
                      eta_grid = if (!is.null(fit$eta_cv)) fit$eta_cv$eta else
                        10^seq(-5, 1, length.out = 12),
                      nfolds = if (!is.null(fit$eta_cv_folds))
                        nrow(fit$eta_cv_folds) else 5L,
                      propensity = NULL, constant_item_effects = NULL,
                      prop_clip = 1e-6, min_valid = 0.90,
                      fit_seed = 1L) {
  X <- as.matrix(X)
  T <- as.numeric(T)
  if (!is.list(Y_list) || !length(Y_list))
    stop("Y_list must be a nonempty list.")
  Y <- lapply(Y_list, as.matrix)
  n <- nrow(X)
  K <- length(Y)
  r <- ncol(Y[[1]])
  if (!is.numeric(X) || ncol(X) < 1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix.")
  if (length(T) != n || anyNA(T) || !all(T %in% 0:1))
    stop("T must contain n binary assignments.")
  if (r < 1L || any(vapply(Y, function(y)
    !is.numeric(y) || ncol(y) != r || any(!is.finite(y)) ||
    !identical(colnames(y), colnames(Y[[1]])), logical(1))))
    stop("Outcome matrices have inconsistent columns or values.")
  ids <- .scoutcrt_ids(X, T, Y, visit_ids)
  if (!is.null(fit$visit_ids) &&
      !identical(ids, lapply(fit$visit_ids, as.integer)))
    stop("visit_ids do not match the SCOUT fit.")
  if (!is.list(fit$final) || length(fit$final) != K ||
      length(fit$alpha) != r || length(fit$eta) != 1L ||
      !is.finite(fit$eta) || fit$eta <= 0)
    stop("Incompatible SCOUT fit.")
  for (k in seq_len(K)) {
    fk <- fit$final[[k]]
    if (!identical(as.integer(fk$ids), ids[[k]]) ||
        nrow(as.matrix(fk$S)) != n ||
        length(fk$CATE) != length(ids[[k]]))
      stop("SCOUT fit rows do not match visit ", k)
  }
  valid_int <- function(z) length(z) == 1L && is.numeric(z) &&
    is.finite(z) && z >= 1 && z == as.integer(z)
  if (!valid_int(B) || !valid_int(cores))
    stop("B and cores must be positive integers.")
  if (!valid_int(nfolds) || nfolds < 2L)
    stop("nfolds must be an integer >= 2.")
  if (!length(eta_grid) || any(!is.finite(eta_grid)) || any(eta_grid <= 0))
    stop("eta_grid must contain positive finite values.")
  if (length(prop_clip) != 1L || !is.finite(prop_clip) ||
      prop_clip <= 0 || prop_clip >= 0.5)
    stop("prop_clip must lie in (0,0.5).")
  if (length(min_valid) != 1L || !is.finite(min_valid) ||
      min_valid <= 0 || min_valid > 1)
    stop("min_valid must lie in (0,1].")
  B <- as.integer(B)
  cores <- as.integer(cores)
  
  labels <- c(paste0("visit", seq_len(K)), "mean_visit_Q")
  J <- length(labels)
  observed <- .scoutcrt_stats(fit)
  S <- do.call(cbind, lapply(fit$final, `[[`, "S"))
  ps <- .scoutcrt_propensity(T, S, prop_clip, propensity)
  delta <- if (is.null(constant_item_effects))
    .scoutcrt_item_effects(X, T, Y, ids, fit$eta) else
      as.matrix(constant_item_effects)
  if (!identical(dim(delta), c(K, r)) || any(!is.finite(delta)))
    stop("constant_item_effects must be a finite K-by-r matrix.")
  
  state <- list(X = X, T = T, Y = Y, ids = ids, delta = delta,
                e = ps$p, n = n, K = K, J = J, seed = seed,
                fit_seed = fit_seed,
                eta = if (isTRUE(reselect_eta)) NULL else fit$eta,
                eta_grid = eta_grid, nfolds = as.integer(nfolds))
  
  if (cores == 1L) {
    draws <- lapply(seq_len(B), .scoutcrt_draw, s = state)
  } else {
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c(
      ".scoutcrt_prep_xt", ".scoutcrt_prep_y", ".scoutcrt_basis",
      ".scoutcrt_nnls", ".scoutcrt_fit", ".scoutcrt_resid",
      ".scoutcrt_R2", ".scoutcrt_refit_pair", ".scoutcrt_stats",
      ".scoutcrt_draw", ".scoutcrt_worker"
    ), envir = environment(SCOUT_CRT))
    parallel::clusterCall(cl, function(st) {
      assign(".scoutcrt_state", st, envir = .GlobalEnv)
      NULL
    }, state)
    draws <- parallel::parLapplyLB(cl, seq_len(B), .scoutcrt_worker,
                                   chunk.size = 5L)
  }
  
  A <- do.call(rbind, lapply(draws, `[[`, "ATE"))
  Q <- do.call(rbind, lapply(draws, `[[`, "Q"))
  colnames(A) <- colnames(Q) <- labels
  p_upper <- function(o, z) (1 + sum(z >= o))/(length(z) + 1)
  
  results <- do.call(rbind, lapply(seq_len(J), function(k) {
    a <- A[is.finite(A[, k]), k]
    q <- Q[is.finite(Q[, k]), k]
    if (min(length(a), length(q)) < min_valid*B)
      stop("Too few successful CRT randomizations at ", labels[k])
    data.frame(
      visit = labels[k],
      n_visit = if (k <= K) length(ids[[k]]) else NA_integer_,
      ATE = observed$ATE[k],
      ATE_p = p_upper(observed$ATE[k], a),
      CATE_Q = observed$Q[k],
      CATE_p = p_upper(observed$Q[k], q)
    )
  }))
  rownames(results) <- NULL
  list(results = results, propensity = ps, constant_item_effects = delta,
       original_fit = fit,
       randomization_distributions = list(ATE = A, Q_homogeneity = Q))
}
