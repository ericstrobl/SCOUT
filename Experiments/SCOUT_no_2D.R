# SCOUT without the final two-dimensional score refinement.
# Supports unequal visit samples via visit_ids; returns the original result fields.
SCOUT_no_2D <- function(X, T, Y, eta = NULL,
                        eta_grid = 10^seq(-5, 1, length.out = 12),
                        nfolds = 5, seed = 1, visit_ids = NULL) {
  X <- as.matrix(X); T <- as.numeric(T); Y <- lapply(Y, as.matrix)
  n <- nrow(X); m <- length(Y); r <- ncol(Y[[1]])
  if (is.null(visit_ids)) visit_ids <- rep(list(seq_len(n)), m)
  if (length(visit_ids) != m || any(vapply(seq_len(m), \(k)
      length(visit_ids[[k]]) != nrow(Y[[k]]), logical(1))))
    stop("visit_ids must match the rows of each outcome matrix")
  visit_ids <- lapply(visit_ids, as.integer)

  res <- function(A, Z) stats::lm.fit(cbind(1, Z), A)$residuals
  R2 <- function(y, Z) {
    Z <- as.matrix(Z); d <- sum((y - mean(y))^2)
    if (!ncol(Z) || d <= .Machine$double.eps) return(0)
    max(0, 1 - sum(res(y, Z)^2) / d)
  }

  prep <- function(X, T, Y, visit_ids) {
    visits <- lapply(seq_along(Y), function(k) {
      ids <- visit_ids[[k]]; Xk <- X[ids, , drop = FALSE]; Tk <- T[ids]
      nk <- length(ids); xm <- colMeans(Xk)
      Xc <- sweep(Xk, 2, xm, "-"); Tc <- Tk - mean(Tk)
      tt <- sum(Tc^2); XX <- crossprod(Xc); xt <- crossprod(Xc, Tc)
      V <- (crossprod(Xc, Y[[k]]) - xt %*% (crossprod(Tc, Y[[k]]) / tt)) / nk
      G <- XX / nk - tcrossprod(drop(xt)) / (2 * nk * tt)
      E <- eigen((G + t(G)) / 2, symmetric = TRUE)
      E$values <- pmax(E$values, 0)
      list(X = Xc, T = Tk, Y = Y[[k]], ids = ids, n = nk,
           p = ncol(X), xm = xm, XX = XX, vT = xt / nk, V = V, E = E)
    })
    list(X = X, T = T, visits = visits, n = nrow(X), m = length(visits))
  }
  solveG <- function(Pk, A, e)
    Pk$E$vectors %*% sweep(crossprod(Pk$E$vectors, A), 1, Pk$E$values + e, "/")
  whiten <- function(Pk, W) {
    C <- crossprod(W, Pk$XX %*% W) / Pk$n
    E <- eigen((C + t(C)) / 2, symmetric = TRUE)
    keep <- E$values > 1e-10 * max(E$values)
    if (!any(keep)) return(matrix(0, Pk$p, 0L))
    W %*% sweep(E$vectors[, keep, drop = FALSE],
                2, sqrt(E$values[keep]), "/")
  }
  adj <- function(A, S) {
    A <- as.matrix(A); A <- sweep(A, 2, colMeans(A), "-")
    if (ncol(S)) A <- A - S %*% (crossprod(S, A) / nrow(S))
    if (ncol(A) == 1L) drop(A) else A
  }

  fit <- function(P, e, effects = TRUE) {
    K <- P$m; pk <- P$visits
    basis <- function(Pk, V) whiten(Pk, solveG(Pk, cbind(Pk$vT, V), e))
    B0 <- lapply(pk, function(Pk) basis(Pk, Pk$V))
    S0_obs <- Map(function(Pk, B) Pk$X %*% B, pk, B0)
    RY <- Map(function(Pk, S) as.matrix(adj(Pk$Y, S)), pk, S0_obs)
    RT <- Map(function(Pk, S) as.numeric(adj(Pk$T, S)), pk, S0_obs)

    nk <- vapply(pk, function(Pk) Pk$n, 0L)
    H <- Reduce(`+`, Map(function(R, n) crossprod(R) / n, RY, nk))
    cvec <- Reduce(`+`, Map(function(R, t, n) drop(crossprod(R, t)) / n,
                            RY, RT, nk))
    tt <- sum(vapply(seq_len(K), function(k) sum(RT[[k]]^2) / nk[k], 0.0))
    fg <- function(a) tt - 2 * sum(a * cvec) + sum(a * drop(H %*% a))
    gr <- function(a) 2 * (drop(H %*% a) - cvec)
    o <- stats::optim(rep(0, r), fg, gr, method = "L-BFGS-B",
                      lower = rep(0, r), control = list(maxit = 500, pgtol = 1e-8))
    a <- pmax(o$par, 0); obj <- fg(a)
    if (sum(a) > 0) a <- a / sum(a)

    if (effects) {
      sd_pool <- sqrt(diag(H) / K)
      marginals <- (cvec / tt) / sd_pool
      names(marginals) <- colnames(pk[[1L]]$Y)
      marginals_by_visit <- do.call(cbind, lapply(seq_len(K), function(k) {
        R <- RY[[k]]; t <- RT[[k]]
        beta <- drop(crossprod(R, t)) / sum(t^2)
        beta / sqrt(colMeans(R^2))
      }))
      rownames(marginals_by_visit) <- colnames(pk[[1L]]$Y)
      colnames(marginals_by_visit) <- paste0("visit_", seq_len(K))
    } else marginals <- marginals_by_visit <- NULL

    final <- lapply(seq_len(K), function(k) {
      Pk <- pk[[k]]; B <- B0[[k]]
      if (!effects) return(list(B = B))
      S_all <- sweep(P$X, 2, Pk$xm, "-") %*% B
      S_obs <- S_all[Pk$ids, , drop = FALSE]
      y <- drop(Pk$Y %*% a); q <- ncol(S_obs)
      Z <- cbind(1, Pk$T, S_obs, S_obs * Pk$T)
      b <- stats::lm.fit(Z, y)$coefficients
      cate <- if (q) drop(b[2] + S_obs %*% tail(b, q)) else rep(b[2], Pk$n)
      cate_all <- if (q) drop(b[2] + S_all %*% tail(b, q)) else rep(b[2], P$n)
      list(B = B, S = S_all, S_visit = S_obs, ids = Pk$ids,
           ATE = mean(cate), CATE = cate, CATE_all = cate_all, coef = b)
    })
    learning_scores <- if (effects) Map(function(Pk, B)
      sweep(P$X, 2, Pk$xm, "-") %*% B, pk, B0) else NULL
    list(alpha = a, objective = obj, learning_scores = learning_scores,
         marginals = marginals, marginals_by_visit = marginals_by_visit,
         final = final)
  }

  if (is.null(eta)) {
    if (any(table(factor(T, levels = 0:1)) < nfolds))
      stop("Each treatment group needs at least nfolds individuals")
    set.seed(seed)
    for (attempt in seq_len(200L)) {
      fold <- integer(n)
      for (g in 0:1) {
        i <- which(T == g)
        fold[i] <- sample(rep(seq_len(nfolds), length.out = length(i)))
      }
      valid <- all(vapply(visit_ids, function(ids)
        all(vapply(seq_len(nfolds), function(f) {
          tr <- ids[fold[ids] != f]; va <- ids[fold[ids] == f]
          length(va) >= 2L && length(unique(T[tr])) == 2L &&
            length(unique(T[va])) == 2L
        }, logical(1))), logical(1)))
      if (valid) break
    }
    if (!valid) stop("Cannot construct stratified folds; reduce nfolds")
    folds <- lapply(seq_len(nfolds), function(f) {
      tr_rows <- lapply(visit_ids, function(ids) which(fold[ids] != f))
      Y_train <- Map(function(y, rows) y[rows, , drop = FALSE], Y, tr_rows)
      ids_train <- Map(function(ids, rows) ids[rows], visit_ids, tr_rows)
      va_rows <- lapply(visit_ids, function(ids) which(fold[ids] == f))
      list(P = prep(X, T, Y_train, ids_train),
           Xval = Map(function(ids, rows) X[ids[rows], , drop = FALSE], visit_ids, va_rows),
           Tval = Map(function(ids, rows) T[ids[rows]], visit_ids, va_rows),
           Yval = Map(function(y, rows) y[rows, , drop = FALSE], Y, va_rows))
    })
    cv <- vapply(eta_grid, function(e) vapply(folds, function(D) {
      M <- fit(D$P, e, effects = FALSE)
      mean(vapply(seq_len(m), function(k) {
        Pk <- D$P$visits[[k]]
        Xval <- sweep(D$Xval[[k]], 2, Pk$xm, "-")
        S <- Xval %*% M$final[[k]]$B
        if (!ncol(S)) return(0)
        y <- drop(D$Yval[[k]] %*% M$alpha); t <- D$Tval[[k]]
        sqrt(R2(t, S) * R2(res(y, t), res(S, t)))
      }, numeric(1)))
    }, numeric(1)), numeric(nfolds))
    mu <- colMeans(cv); eta <- eta_grid[which.max(mu)]
    eta_cv <- data.frame(eta = eta_grid, mean_score = mu,
                         se_score = apply(cv, 2, sd) / sqrt(nfolds))
  } else cv <- eta_cv <- NULL

  M <- fit(prep(X, T, Y, visit_ids), eta, effects = TRUE)
  list(alpha = M$alpha, active = which(M$alpha > 1e-8),
       objective = M$objective, marginals = M$marginals,
       marginals_by_visit = M$marginals_by_visit,
       eta = eta, eta_cv = eta_cv, eta_cv_folds = cv,
       learning_scores = M$learning_scores, final = M$final,
       visit_ids = visit_ids)
}
