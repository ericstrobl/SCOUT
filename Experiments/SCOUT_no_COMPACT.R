#' SCOUT ablation without COMPACT adjustment
#'
#' `SCOUT_no_COMPACT()` replaces COMPACT's low-dimensional proxy adjustment
#' with ridge regression on the full baseline covariates. A single nonnegative
#' composite outcome is learned across visits, followed by a ridge treatment-
#' effect model with baseline covariates and treatment-by-covariate interactions.
#' This isolates the contribution of COMPACT within the SCOUT framework.
#'
#' Visits may have different observed individuals. `X` and `T` contain all `n`
#' individuals; rows of `Y[[k]]` correspond to `visit_ids[[k]]` in that order.
#' Treatment is assigned once per individual and shared across visits.
#'
#' @param X An `n` by `p` numeric matrix of finite baseline covariates.
#' @param T Binary numeric vector of length `n` (1 = treatment, 0 = comparator).
#' @param Y List of `m` numeric outcome matrices, with `n_k` individuals and
#'   the same `r` outcome components in each matrix.
#' @param visit_ids Optional list of `m` distinct baseline row-index vectors,
#'   matching the rows of `Y[[k]]` in order. If `NULL`, all visits contain all
#'   `n` individuals in baseline order.
#' @param eta Ridge penalty used to adjust treatment and outcomes during
#'   outcome learning. If `NULL`, select it by subject-level cross-validation
#'   of the held-out ridge-adjusted residual loss for the learned composite.
#' @param eta_grid Positive candidate ridge penalties. Default:
#'   `10^seq(-5, 1, length.out = 12)`.
#' @param nfolds Number of cross-validation folds for both `eta` selection and
#'   each final ridge treatment-effect model. Must be at least 3. Folds for
#'   `eta` selection are assigned by subject across visits. Default: `5`.
#' @param seed Integer random seed for cross-validation fold construction.
#'   Default: `1`.
#'
#' @details
#' At visit `k`, treatment and each outcome component are ridge-adjusted for
#' the original `p` baseline covariates. A common nonnegative outcome-weight
#' vector minimizes the sum of residual squared errors weighted by `1/n_k`,
#' giving every visit equal weight. The NNLS solution is normalized to sum to
#' one if nonzero; `objective` refers to the unnormalized NNLS optimum.
#'
#' Unlike `SCOUT()`, this ablation does not learn COMPACT adjustment scores.
#' Each final visit-specific model uses ridge regression of the learned outcome
#' on treatment, centered baseline covariates, and treatment interactions.
#' Treatment is unpenalized; covariate and interaction coefficients are ridge-
#' penalized. `lambda.min` is chosen by cross-validation at each visit.
#'
#' @return A list with components:
#' \describe{
#'   \item{alpha}{Normalized nonnegative outcome weights (length `r`).}
#'   \item{active}{Indices of weights greater than `1e-8`.}
#'   \item{objective}{Visit-normalized NNLS objective before weight scaling.}
#'   \item{marginals}{Pooled standardized ridge-adjusted marginal associations.}
#'   \item{marginals_by_visit}{Visit-specific adjusted associations (`r` by `m`).}
#'   \item{eta}{Selected or supplied ridge adjustment penalty.}
#'   \item{lambda}{Visit-specific penalties selected for final ridge models.}
#'   \item{eta_cv}{Cross-validation summary or `NULL` if `eta` was supplied.}
#'   \item{eta_cv_folds}{Fold-by-parameter cross-validation scores or `NULL`.}
#'   \item{learning_scores}{`NULL`; no COMPACT score space is used.}
#'   \item{final}{List of visit-specific results: `ids`, `ATE`, `CATE`
#'     (length `n_k`), `CATE_all` (length `n`), `coef`, and `lambda`.
#'     `S` and `S_visit` are `NULL` because the full covariates are used.}
#'   \item{visit_ids}{Validated observation-index vectors for each visit.}
#' }
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' X <- matrix(rnorm(n * 8), n, 8)
#' T <- rbinom(n, 1, 0.5)
#' ids <- list(seq_len(n), sort(sample(seq_len(n), 55)))
#' Y <- lapply(ids, function(i) {
#'   cbind(X[i, 1] + 0.5 * T[i] + rnorm(length(i)),
#'         X[i, 2] + T[i] + rnorm(length(i)))
#' })
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   fit <- SCOUT_no_COMPACT(X, T, Y, eta = 0.1, visit_ids = ids)
#'   fit$alpha
#'   fit$final[[2]]$ATE
#' }
#'
#' @seealso [SCOUT()], [SCOUT_no_2D()]
#' @md
#' @export
SCOUT_no_COMPACT <- function(X, T, Y, eta = NULL,
                       eta_grid = 10^seq(-5, 1, length.out = 12),
                       nfolds = 5, seed = 1, visit_ids = NULL) {

  X <- as.matrix(X)
  T <- as.numeric(T)
  if (!is.numeric(X) || any(!is.finite(X)) || ncol(X) < 1L)
    stop("X must be a finite numeric matrix with at least one column.")
  n <- nrow(X)
  if (length(T) != n || anyNA(T) || !all(T %in% 0:1))
    stop("T must contain one binary treatment assignment per row of X.")
  if (!is.list(Y) || !length(Y))
    stop("Y must be a nonempty list of outcome matrices.")
  Y <- lapply(Y, as.matrix)
  if (any(vapply(Y, function(y) !identical(colnames(y), colnames(Y[[1L]])), logical(1))))
    stop("All Y matrices must have identical outcome column names and ordering.")
  m <- length(Y)
  r <- ncol(Y[[1L]])
  if (!is.numeric(Y[[1L]]) || r < 1L ||
      any(vapply(Y, function(y) !is.numeric(y) || ncol(y) != r ||
                   any(!is.finite(y)), logical(1))))
    stop("All Y matrices must contain finite numeric outcomes in the same column order.")

  # Each Y[[k]] row corresponds to visit_ids[[k]] in the full baseline X/T.
  if (is.null(visit_ids)) {
    if (any(vapply(Y, nrow, 0L) != n))
      stop("Supply visit_ids when a visit contains fewer than n individuals.")
    visit_ids <- rep(list(seq_len(n)), m)
  }
  if (!is.list(visit_ids) || length(visit_ids) != m)
    stop("visit_ids must be a list with one index vector per Y matrix.")
  for (k in seq_len(m)) {
    ids <- visit_ids[[k]]
    if (!is.numeric(ids) || anyNA(ids) || any(!is.finite(ids)) ||
        any(ids != as.integer(ids)) || any(ids < 1L | ids > n) ||
        anyDuplicated(ids) || length(ids) != nrow(Y[[k]]))
      stop("Invalid visit_ids[[", k, "]]: provide distinct baseline row indices matching Y[[", k, "]] rows.")
    if (length(ids) < 3L || length(unique(T[ids])) < 2L)
      stop("Visit ", k, " needs at least three individuals and both treatment groups.")
    visit_ids[[k]] <- as.integer(ids)
  }
  if (!is.null(eta) && (length(eta) != 1L || !is.finite(eta) || eta <= 0))
    stop("eta must be NULL or a positive finite scalar.")
  if (length(eta_grid) < 1L || any(!is.finite(eta_grid) | eta_grid <= 0))
    stop("eta_grid must contain positive finite values.")
  if (length(nfolds) != 1L || !is.finite(nfolds) ||
      nfolds != as.integer(nfolds) || nfolds < 3L)
    stop("nfolds must be an integer greater than or equal to three.")

  # Precompute visit-specific ridge quantities using only observed subjects.
  prep <- function(X, T, Y, visit_ids) {
    visits <- lapply(seq_along(Y), function(k) {
      ids <- visit_ids[[k]]
      Xk <- X[ids, , drop = FALSE]
      Tk <- T[ids]
      nk <- length(ids)
      xm <- colMeans(Xk)
      Xc <- sweep(Xk, 2, xm, "-")
      Tc <- Tk - mean(Tk)
      XX <- crossprod(Xc) / nk
      E <- eigen((XX + t(XX)) / 2, symmetric = TRUE)
      E$values <- pmax(E$values, 0)
      list(X = Xc, T = Tk, Y = Y[[k]], ids = ids, n = nk,
           p = ncol(X), xm = xm, vT = crossprod(Xc, Tc) / nk,
           U = crossprod(Xc, Y[[k]]) / nk, E = E)
    })
    list(X = X, T = T, visits = visits, n = nrow(X), m = length(visits))
  }

  # Apply (X_k'X_k/n_k + eta I)^(-1) without forming a matrix inverse.
  solveR <- function(Pk, A, e)
    Pk$E$vectors %*%
      sweep(crossprod(Pk$E$vectors, A), 1, Pk$E$values + e, "/")

  fit <- function(P, e, effects = TRUE) {
    K <- P$m
    pk <- P$visits

    # Stage 1: ridge-adjust treatment and each outcome against all p covariates.
    bT <- lapply(pk, function(Pk) solveR(Pk, Pk$vT, e))
    bY <- lapply(pk, function(Pk) solveR(Pk, Pk$U, e))
    RT <- Map(function(Pk, b) as.numeric(Pk$T - mean(Pk$T) - Pk$X %*% b),
              pk, bT)
    RY <- Map(function(Pk, b) sweep(Pk$Y, 2, colMeans(Pk$Y), "-") - Pk$X %*% b,
              pk, bY)

    # Stage 2: same shared-outcome NNLS as SCOUT, with visit weights 1/n_k.
    nk <- vapply(pk, function(Pk) Pk$n, 0L)
    H <- Reduce(`+`, Map(function(R, n) crossprod(R) / n, RY, nk))
    cvec <- Reduce(`+`, Map(function(R, t, n) drop(crossprod(R, t)) / n,
                            RY, RT, nk))
    tt <- sum(vapply(seq_len(K), function(k) sum(RT[[k]]^2) / nk[k], 0.0))
    if (tt <= .Machine$double.eps)
      stop("No residual treatment variation remains after ridge adjustment.")
    fg <- function(a) tt - 2 * sum(a * cvec) + sum(a * drop(H %*% a))
    gr <- function(a) 2 * (drop(H %*% a) - cvec)
    o <- stats::optim(rep(0, r), fg, gr, method = "L-BFGS-B",
                      lower = rep(0, r),
                      control = list(maxit = 500, pgtol = 1e-8))
    a <- pmax(o$par, 0)
    obj <- fg(a)  # Record the NNLS optimum before sum-to-one normalization.
    if (sum(a) > 0) a <- a / sum(a)

    # For eta selection: predict the ridge-adjusted T and learned outcome
    # on new subjects, using the training visit's centering constants.
    if (!effects) {
      return(list(alpha = a, objective = obj, bT = bT,
                  bY = lapply(pk, function(Pk) solveR(Pk, Pk$U %*% a, e)),
                  tm = vapply(pk, function(Pk) mean(Pk$T), 0.0),
                  ym = vapply(pk, function(Pk) mean(drop(Pk$Y %*% a)), 0.0)))
    }

    # Descriptive ridge-adjusted marginal treatment/outcome associations.
    sd_pool <- sqrt(diag(H) / K)
    marginals <- (cvec / tt) / sd_pool
    names(marginals) <- colnames(pk[[1L]]$Y)
    marginals_by_visit <- do.call(cbind, lapply(seq_len(K), function(k) {
      R <- RY[[k]]
      t <- RT[[k]]
      beta <- drop(crossprod(R, t)) / sum(t^2)
      beta / sqrt(colMeans(R^2))
    }))
    rownames(marginals_by_visit) <- colnames(pk[[1L]]$Y)
    colnames(marginals_by_visit) <- paste0("visit_", seq_len(K))

    # Stage 3: ridge CATE estimation using X and T:X, not COMPACT scores.
    # Fit each visit only on its observed n_k subjects. Use the same fitted
    # model to predict CATE for all baseline individuals if requested.
    final <- lapply(seq_len(K), function(k) {
      Pk <- pk[[k]]
      nk <- Pk$n
      y <- drop(Pk$Y %*% a)
      Z <- cbind(Pk$T, Pk$X, Pk$X * Pk$T)

      ng <- table(factor(Pk$T, levels = 0:1))
      Kfold <- min(nfolds, ng)
      if (Kfold < 3L)
        stop("At least 3 subjects per treatment group are needed at visit ", k,
             " for ridge-model cross-validation.")
      set.seed(seed + 100000L + k)
      foldid <- integer(nk)
      for (g in 0:1) {
        i <- which(Pk$T == g)
        foldid[i] <- sample(rep(seq_len(Kfold), length.out = length(i)))
      }
      f <- glmnet::cv.glmnet(Z, y, alpha = 0, foldid = foldid,
                             type.measure = "mse", standardize = TRUE,
                             penalty.factor = c(0, rep(1, 2 * Pk$p)))
      b <- as.numeric(stats::coef(f, s = "lambda.min"))
      bx <- b[seq.int(Pk$p + 3L, 2L * Pk$p + 2L)]
      cate <- drop(b[2L] + Pk$X %*% bx)
      X_all <- sweep(P$X, 2, Pk$xm, "-")
      cate_all <- drop(b[2L] + X_all %*% bx)

      list(B = NULL, S = NULL, S_visit = NULL, ids = Pk$ids,
           ATE = mean(cate), CATE = cate, CATE_all = cate_all,
           coef = b, lambda = f$lambda.min)
    })

    list(alpha = a, objective = obj, learning_scores = NULL,
         marginals = marginals, marginals_by_visit = marginals_by_visit,
         final = final)
  }

  # Choose eta by subject-level cross-validation. A patient appearing in more
  # than one visit has one fold assignment, preventing leakage across visits.
  if (is.null(eta)) {
    if (any(table(factor(T, levels = 0:1)) < nfolds))
      stop("Each treatment group needs at least nfolds baseline individuals.")
    if (any(vapply(visit_ids, function(ids)
          any(table(factor(T[ids], levels = 0:1)) < nfolds), logical(1))))
      stop("Each treatment group at each visit needs at least nfolds individuals for CV.")

    set.seed(seed)
    # Re-draw global stratified folds if necessary to ensure both treatments
    # are represented among training and validation subjects at every visit.
    for (attempt in seq_len(200L)) {
      fold <- integer(n)
      for (g in 0:1) {
        i <- which(T == g)
        fold[i] <- sample(rep(seq_len(nfolds), length.out = length(i)))
      }
      valid <- all(vapply(visit_ids, function(ids) {
        all(vapply(seq_len(nfolds), function(f) {
          tr <- ids[fold[ids] != f]
          va <- ids[fold[ids] == f]
          length(va) >= 2L && length(unique(T[tr])) == 2L &&
            length(unique(T[va])) == 2L
        }, logical(1)))
      }, logical(1)))
      if (valid) break
    }
    if (!valid)
      stop("Cannot form subject-level CV folds with both treatments at every visit; ",
           "reduce nfolds or increase per-visit sample sizes.")

    folds <- lapply(seq_len(nfolds), function(f) {
      # Restrict each outcome matrix and index list to training subjects.
      tr_rows <- lapply(visit_ids, function(ids) which(fold[ids] != f))
      Y_train <- Map(function(y, rows) y[rows, , drop = FALSE], Y, tr_rows)
      ids_train <- Map(function(ids, rows) ids[rows], visit_ids, tr_rows)

      # Validation is similarly restricted *within each visit*.
      va_rows <- lapply(visit_ids, function(ids) which(fold[ids] == f))
      list(P = prep(X, T, Y_train, ids_train),
           Xval = Map(function(ids, rows) X[ids[rows], , drop = FALSE],
                      visit_ids, va_rows),
           Tval = Map(function(ids, rows) T[ids[rows]], visit_ids, va_rows),
           Yval = Map(function(y, rows) y[rows, , drop = FALSE], Y, va_rows))
    })

    # Unlike COMPACT, ridge has no learned low-dimensional score space.
    # Select eta using the original ablation's held-out residual loss.
    cv <- vapply(eta_grid, function(e) vapply(folds, function(D) {
      M <- fit(D$P, e, effects = FALSE)
      mean(vapply(seq_len(m), function(k) {
        Pk <- D$P$visits[[k]]
        Xval <- sweep(D$Xval[[k]], 2, Pk$xm, "-")
        rt <- drop(D$Tval[[k]] - M$tm[k] - Xval %*% M$bT[[k]])
        y <- drop(D$Yval[[k]] %*% M$alpha)
        ry <- drop(y - M$ym[k] - Xval %*% M$bY[[k]])
        -mean((rt - ry)^2)
      }, numeric(1)))
    }, numeric(1)), numeric(nfolds))

    mu <- colMeans(cv)
    eta <- eta_grid[which.max(mu)]
    eta_cv <- data.frame(eta = eta_grid, mean_score = mu,
                         se_score = apply(cv, 2, sd) / sqrt(nfolds))
  } else {
    cv <- eta_cv <- NULL
  }

  # Fit on all available outcomes at the selected eta.
  M <- fit(prep(X, T, Y, visit_ids), eta, effects = TRUE)
  list(alpha = M$alpha,
       active = which(M$alpha > 1e-8),
       objective = M$objective,
       marginals = M$marginals,
       marginals_by_visit = M$marginals_by_visit,
       eta = eta,
       lambda = vapply(M$final, function(x) x$lambda, 0.0),
       eta_cv = eta_cv,
       eta_cv_folds = cv,
       learning_scores = M$learning_scores,
       final = M$final,
       visit_ids = visit_ids)
}
