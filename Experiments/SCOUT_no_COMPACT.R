#' SCOUT without COMPACT adjustment
#'
#' Ridge-adjusts treatment and outcomes using all baseline covariates, learns
#' nonnegative outcome weights, and fits ridge treatment-interaction models.
#' For complete visits, follows the original SCOUT_NNLS_ridge computation;
#' outcome weights and effects are normalized only after the final fits.
#' @export
SCOUT_no_COMPACT <- function(X, T, Y, eta = NULL,
                             eta_grid = 10^seq(-5, 1, length.out = 12),
                             nfolds = 5, seed = 1, visit_ids = NULL) {
  
  X <- as.matrix(X); T <- as.numeric(T); Y <- lapply(Y, as.matrix)
  n <- nrow(X); p <- ncol(X); m <- length(Y); r <- ncol(Y[[1]])
  if (is.null(visit_ids)) visit_ids <- rep(list(seq_len(n)), m)
  if (length(visit_ids) != m || any(vapply(seq_len(m), \(k)
                                           length(visit_ids[[k]]) != nrow(Y[[k]]), logical(1))))
    stop("visit_ids must match the rows of each outcome matrix")
  
  # Reuse one decomposition when the same subjects appear at every visit.
  prep <- function(Ys, ids) {
    one <- function(i, y) {
      x <- X[i, , drop = FALSE]; t <- T[i]; nk <- length(i)
      xm <- colMeans(x); x <- sweep(x, 2, xm, "-")
      tc <- t - mean(t); xx <- crossprod(x); xt <- crossprod(x, tc)
      u <- crossprod(x, y) / nk
      E <- eigen(xx / nk, symmetric = TRUE)
      E$values <- pmax(E$values, 0)
      list(X = x, T = t, Y = y, n = nk, p = p, xm = xm,
           vT = xt / nk, U = u, E = E, ids = i)
    }
    shared <- all(vapply(ids, identical, logical(1), ids[[1]]))
    if (shared) {
      first <- one(ids[[1]], Ys[[1]])
      visits <- lapply(seq_len(m), \(k) {
        if (k == 1L) return(first)
        q <- first; q$Y <- Ys[[k]]
        q$U <- crossprod(q$X, q$Y) / q$n
        q
      })
    } else {
      visits <- Map(one, ids, Ys)
    }
    list(visits = visits, shared = shared)
  }
  
  solveR <- function(P, A, e)
    P$E$vectors %*%
    sweep(crossprod(P$E$vectors, A), 1, P$E$values + e, "/")
  
  fit <- function(P, e, effects = TRUE) {
    pk <- P$visits
    ridge <- function(Q, A, V) {
      A <- as.matrix(A); A <- sweep(A, 2, colMeans(A), "-")
      Z <- A - Q$X %*% solveR(Q, V, e)
      if (ncol(Z) == 1L) drop(Z) else Z
    }
    RT <- if (P$shared) {
      rep(list(ridge(pk[[1]], pk[[1]]$T, pk[[1]]$vT)), m)
    } else lapply(pk, \(q) ridge(q, q$T, q$vT))
    RY <- lapply(pk, \(q) ridge(q, q$Y, q$U))
    
    # Original unscaled NNLS (never normalize before eta CV or glmnet).
    H <- Reduce(`+`, lapply(RY, crossprod))
    c <- Reduce(`+`, Map(\(y, t) drop(crossprod(y, t)), RY, RT))
    tt <- if (P$shared) m * sum(RT[[1]]^2) else
      sum(vapply(RT, \(t) sum(t^2), 0.0))
    fg <- function(a) tt - 2 * sum(a * c) + sum(a * drop(H %*% a))
    gr <- function(a) 2 * (drop(H %*% a) - c)
    o <- optim(rep(0, r), fg, gr, method = "L-BFGS-B",
               lower = rep(0, r), control = list(maxit = 500, pgtol = 1e-8))
    a <- pmax(o$par, 0)
    
    bT <- lapply(pk, \(q) solveR(q, q$vT, e))
    bY <- lapply(pk, \(q) solveR(q, q$U %*% a, e))
    ym <- vapply(pk, \(q) mean(drop(q$Y %*% a)), 0.0)
    if (!effects) return(list(alpha = a, bT = bT, bY = bY,
                              tm = vapply(pk, \(q) mean(q$T), 0.0), ym = ym))
    
    final <- lapply(seq_len(m), \(k) {
      q <- pk[[k]]
      y <- drop(q$Y %*% a)
      Z <- cbind(q$T, q$X, q$X * q$T)
      K <- min(nfolds, min(table(q$T)))
      set.seed(seed + 100000 + k)
      id <- integer(q$n)
      for (g in unique(q$T)) {
        i <- which(q$T == g)
        id[i] <- sample(rep(seq_len(K), length.out = length(i)))
      }
      f <- glmnet::cv.glmnet(Z, y, alpha = 0, foldid = id,
                             type.measure = "mse",
                             penalty.factor = c(0, rep(1, 2 * p)),
                             standardize = TRUE)
      b <- as.numeric(coef(f, s = "lambda.min"))
      bx <- b[(p + 3):(2 * p + 2)]
      cate <- drop(b[2] + q$X %*% bx)
      cate_all <- if (P$shared) cate else
        drop(b[2] + sweep(X, 2, q$xm, "-") %*% bx)
      list(ids = q$ids, S = NULL, S_visit = NULL,
           ATE = mean(cate), CATE = cate, CATE_all = cate_all,
           coef = b, lambda = f$lambda.min)
    })
    list(alpha = a, objective = fg(a), final = final)
  }
  
  if (is.null(eta)) {
    if (any(table(T) < nfolds)) stop("Each treatment group needs >= nfolds")
    set.seed(seed)
    fold <- integer(n)
    for (g in unique(T)) {
      i <- which(T == g)
      fold[i] <- sample(rep(seq_len(nfolds), length.out = length(i)))
    }
    F <- lapply(seq_len(nfolds), \(k) {
      train <- lapply(visit_ids, \(i) fold[i] != k)
      train_ids <- Map(\(i, s) i[s], visit_ids, train)
      test_ids <- Map(\(i, s) i[!s], visit_ids, train)
      train_Y <- Map(\(y, s) y[s, , drop = FALSE], Y, train)
      test_Y <- Map(\(y, s) y[!s, , drop = FALSE], Y, train)
      list(P = prep(train_Y, train_ids), ids = test_ids, Y = test_Y)
    })
    cv <- vapply(eta_grid, \(e) vapply(F, \(D) {
      M <- fit(D$P, e, FALSE)
      mean(vapply(seq_len(m), \(k) {
        q <- D$P$visits[[k]]
        xv <- sweep(X[D$ids[[k]], , drop = FALSE], 2, q$xm, "-")
        rt <- drop(T[D$ids[[k]]] - M$tm[k] - xv %*% M$bT[[k]])
        y <- drop(D$Y[[k]] %*% M$alpha)
        ry <- drop(y - M$ym[k] - xv %*% M$bY[[k]])
        -mean((rt - ry)^2)
      }, numeric(1)))
    }, numeric(1)), numeric(nfolds))
    mu <- colMeans(cv)
    eta <- eta_grid[which.max(mu)]
    eta_cv <- data.frame(eta = eta_grid, mean_score = mu,
                         se_score = apply(cv, 2, sd) / sqrt(nfolds))
  } else cv <- eta_cv <- NULL
  
  # Fit the ridge ablation, then normalize.
  M <- fit(prep(Y, visit_ids), eta, TRUE)
  s <- sum(abs(M$alpha))
  if (!is.finite(s) || s == 0) stop("NNLS returned zero outcome weights")
  alpha_raw <- M$alpha
  M$alpha <- alpha_raw / s
  M$final <- lapply(M$final, \(f) {
    f$ATE <- f$ATE / s
    f$CATE <- f$CATE / s
    f$CATE_all <- f$CATE_all / s
    f
  })
  c(M, list(alpha_raw = alpha_raw, normalization_factor = s,
            active = which(alpha_raw > 1e-8), eta = eta,
            lambda = vapply(M$final, \(f) f$lambda, 0.0),
            eta_cv = eta_cv, eta_cv_folds = cv, visit_ids = visit_ids))
}
