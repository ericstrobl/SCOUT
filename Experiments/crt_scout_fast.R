# -----------------------------------------------------------------------------
# Fast paired SCOUT refits: equivalent complete-visit NNLS/COMPACT calculation
# without repeatedly constructing n-by-r residual matrices during eta CV.
# -----------------------------------------------------------------------------
crt_scout_multi <- function(X, T, Ysets, eta = NULL,
                            eta_grid = 10^seq(-5, 1, length.out = 12),
                            nfolds = 5, seed = 1,
                            return_learning_scores = FALSE,
                            estimate_item_effects = FALSE) {
  X <- as.matrix(X); T <- as.numeric(T)
  Ysets <- lapply(Ysets, \(Y) lapply(Y, as.matrix))
  nd <- length(Ysets); K <- length(Ysets[[1]]); r <- ncol(Ysets[[1]][[1]])
  if (!nd || any(vapply(Ysets, length, 0L) != K)) stop("Inconsistent outcome sets")
  if (length(T) != nrow(X) || anyNA(T) || !all(T %in% 0:1)) stop("Invalid T")
  if (any(vapply(unlist(Ysets, recursive = FALSE), \(y)
                 nrow(y) != length(T) || ncol(y) != r || any(!is.finite(y)), TRUE)))
    stop("Outcome dimensions or values invalid")

  # Treatment-dependent matrices are reused across both outcome sets.
  prep_XT <- function(X, T) {
    n <- nrow(X); p <- ncol(X); xm <- colMeans(X)
    X <- sweep(X, 2, xm, "-"); tc <- T - mean(T); tt <- sum(tc^2)
    if (tt <= 0) stop("Treatment has no variation")
    XX <- crossprod(X); xt <- crossprod(X, tc)
    G <- XX/n - tcrossprod(drop(xt))/(2*n*tt)
    E <- eigen((G + t(G))/2, symmetric = TRUE)
    list(X = X, T = T, tc = tc, tt = tt, n = n, p = p, xm = xm,
         U = E$vectors, raw_ev = E$values, ev = pmax(E$values, 0),
         uT = crossprod(E$vectors, xt/n))
  }

  prep_Y <- function(P, Y) {
    YC <- lapply(Y, \(y) sweep(y, 2, colMeans(y), "-"))
    TY <- lapply(YC, \(y) crossprod(P$tc, y))
    ZY <- lapply(YC, \(y) crossprod(P$U, crossprod(P$X, y)/P$n))
    V <- Map(\(zy, ty) zy - P$uT %*% (ty/P$tt), ZY, TY)
    list(Y = Y, V = V, ZY = ZY, TY = TY, YY = lapply(YC, crossprod))
  }

  # Whiten adjustment directions in the eigenspace of G, using only small
  # matrices. Identity: U'XXU/n = diag(raw_ev) + n*uT*uT'/(2*tt).
  score_basis <- function(P, Q, e, need_B = FALSE) {
    D <- sweep(Q, 1, P$ev + e, "/")
    C <- crossprod(D, sweep(D, 1, P$raw_ev, "*"))
    a <- drop(crossprod(D, P$uT))
    C <- C + (P$n/(2*P$tt))*tcrossprod(a)
    E <- eigen((C + t(C))/2, symmetric = TRUE)
    keep <- E$values > 1e-10*max(E$values)
    if (!any(keep)) {
      R <- matrix(0, P$p, 0L); gram <- matrix(0, 0L, 0L)
    } else {
      W <- sweep(E$vectors[, keep, drop = FALSE], 2, sqrt(E$values[keep]), "/")
      R <- D %*% W
      gram <- crossprod(W, C %*% W)
    }
    list(R = R, gram = gram, B = if (need_B) P$U %*% R else NULL)
  }

  # Estimate constant effects for each outcome component using its own
  # treatment + outcome-specific two-score COMPACT adjustment.
  item_effects <- function(P, O, e) {
    D <- cbind(P$uT, do.call(cbind, O$V))
    S <- P$X %*% (P$U %*% sweep(D, 1, P$ev + e, "/"))
    delta <- matrix(NA_real_, K, r)
    for (k in seq_len(K)) for (j in seq_len(r)) {
      h <- 1L + (k - 1L)*r + j
      delta[k, j] <- lm.fit(cbind(1, P$T, S[, 1L], S[, h]), O$Y[[k]][, j])$coefficients[2L]
    }
    if (any(!is.finite(delta))) stop("Nonfinite outcome-specific COMPACT effects")
    delta
  }

  nnls <- function(H, c) {
    f <- function(a) sum(a*drop(H %*% a)) - 2*sum(a*c)
    g <- function(a) 2*(drop(H %*% a) - c)
    pmax(optim(rep(0, length(c)), f, g, method = "L-BFGS-B",
               lower = rep(0, length(c)),
               control = list(maxit = 500, pgtol = 1e-8))$par, 0)
  }

  fit_one <- function(P, O, e, effects = TRUE, learning = FALSE) {
    # Sufficient statistics for the visit-weighted NNLS outcome objective.
    H <- matrix(0, r, r); c <- numeric(r); tt <- 0
    for (k in seq_len(K)) {
      z <- score_basis(P, cbind(P$uT, O$V[[k]]), e)
      if (!ncol(z$R)) {
        H <- H + O$YY[[k]]
        c <- c + drop(O$TY[[k]])
        tt <- tt + P$tt
      } else {
        sy <- crossprod(z$R, O$ZY[[k]])
        st <- drop(crossprod(z$R, P$uT))
        G <- z$gram
        H <- H + O$YY[[k]] - 2*P$n*crossprod(sy) +
          P$n*crossprod(sy, G %*% sy)
        c <- c + drop(O$TY[[k]]) - 2*P$n*drop(crossprod(sy, st)) +
          P$n*drop(crossprod(sy, G %*% st))
        tt <- tt + P$tt - 2*P$n*sum(st^2) +
          P$n*drop(crossprod(st, G %*% st))
      }
    }
    H <- (H + t(H))/2
    alpha_raw <- nnls(H/P$n, c/P$n)
    alpha <- alpha_raw
    if (sum(alpha) > 0) alpha <- alpha/sum(alpha)

    # Final visit-specific, at-most-2D COMPACT scores and treatment effects.
    final <- lapply(seq_len(K), \(k) {
      z <- score_basis(P, cbind(P$uT, O$V[[k]] %*% alpha), e, need_B = TRUE)
      if (!effects) return(list(B = z$B))
      S <- P$X %*% z$B; y <- drop(O$Y[[k]] %*% alpha); q <- ncol(S)
      b <- lm.fit(cbind(1, P$T, S, S*P$T), y)$coefficients
      cate <- if (q) drop(b[2] + S %*% tail(b, q)) else rep(b[2], P$n)
      list(B = z$B, S = S, ATE = mean(cate), CATE = cate, coef = b)
    })
    S0 <- if (effects && learning) lapply(seq_len(K), \(k) {
      z <- score_basis(P, cbind(P$uT, O$V[[k]]), e, need_B = TRUE)
      P$X %*% z$B
    }) else NULL

    # Store objective at the raw NNLS optimum, before alpha normalization.
    objective <- if (effects) tt - 2*sum(alpha_raw*c) +
      sum(alpha_raw*drop(H %*% alpha_raw)) else NA_real_
    list(alpha = alpha, objective = objective, learning_scores = S0,
         final = final, H = H/P$n, c = c/P$n, n = P$n)
  }

  R2 <- function(y, Z) {
    Z <- as.matrix(Z); d <- sum((y - mean(y))^2)
    if (!ncol(Z) || d <= .Machine$double.eps) return(0)
    max(0, 1 - sum(lm.fit(cbind(1, Z), y)$residuals^2)/d)
  }
  residual <- function(A, Z) lm.fit(cbind(1, Z), A)$residuals

  if (is.null(eta)) {
    if (any(table(T) < nfolds)) stop("Each treatment group needs >= nfolds")
    set.seed(seed); ids <- integer(length(T))
    for (g in unique(T)) {
      ix <- which(T == g)
      ids[ix] <- sample(rep(seq_len(nfolds), length.out = length(ix)))
    }
    folds <- lapply(seq_len(nfolds), \(k) {
      tr <- ids != k; P <- prep_XT(X[tr, , drop = FALSE], T[tr])
      list(P = P, Xtest = sweep(X[!tr, , drop = FALSE], 2, P$xm, "-"),
           Ttest = T[!tr],
           train = lapply(Ysets, \(Y) prep_Y(P, lapply(Y, \(y) y[tr, , drop = FALSE]))),
           test = lapply(Ysets, \(Y) lapply(Y, \(y) y[!tr, , drop = FALSE])))
    })
    # Independently choose eta for every outcome set using the same fold splits.
    CV <- lapply(seq_len(nd), \(i) matrix(NA_real_, nfolds, length(eta_grid)))
    for (k in seq_len(nfolds)) for (j in seq_along(eta_grid)) {
      D <- folds[[k]]
      for (i in seq_len(nd)) {
        f <- fit_one(D$P, D$train[[i]], eta_grid[j], FALSE)
        CV[[i]][k, j] <- mean(vapply(seq_len(K), \(v) {
          S <- D$Xtest %*% f$final[[v]]$B
          y <- drop(D$test[[i]][[v]] %*% f$alpha)
          sqrt(R2(D$Ttest, S)*R2(residual(y, D$Ttest), residual(S, D$Ttest)))
        }, numeric(1)))
      }
    }
    means <- lapply(CV, colMeans)
    etas <- vapply(means, \(z) eta_grid[which.max(z)], numeric(1))
    eta_cv <- Map(\(mu, cv) data.frame(eta = eta_grid, mean_score = mu,
                                      se_score = apply(cv, 2, sd)/sqrt(nfolds)), means, CV)
  } else {
    etas <- rep(eta, nd); CV <- eta_cv <- rep(list(NULL), nd)
  }

  # Share the full-sample treatment eigensystem between the two null fits.
  P <- prep_XT(X, T)
  lapply(seq_len(nd), \(i) {
    O <- prep_Y(P, Ysets[[i]])
    f <- fit_one(P, O, etas[i], TRUE, return_learning_scores)
    list(alpha = f$alpha, active = which(f$alpha > 1e-8),
         objective = f$objective, eta = etas[i], eta_cv = eta_cv[[i]],
         eta_cv_folds = CV[[i]], learning_scores = f$learning_scores,
         final = f$final, H = f$H, c = f$c, n = f$n,
         constant_item_effects = if (estimate_item_effects)
           item_effects(P, O, etas[i]) else NULL)
  })
}

crt_scout <- function(X, T, Y, eta = NULL, seed = 1, estimate_item_effects = FALSE)
  crt_scout_multi(X, T, list(Y), eta = eta, seed = seed,
                  return_learning_scores = TRUE,
                  estimate_item_effects = estimate_item_effects)[[1]]

