#' SCOUT ablation without the final two-dimensional refinement
#'
#' `SCOUT_no_2D()` learns a shared, nonnegative composite of longitudinal
#' outcome components while adjusting for high-dimensional, noisy baseline
#' covariates through COMPACT scores. Unlike `SCOUT()`, it retains the initial
#' (at most `r + 1`)-dimensional adjustment scores for final ATE and CATE
#' estimation rather than refining them to two outcome-specific scores.
#' All other steps use the same procedure as `SCOUT()`.
#'
#' Different visits may contain different sets of observed individuals. The
#' baseline matrix `X` and binary treatment vector `T` contain all `n`
#' individuals, while `Y[[k]]` contains only those observed at visit `k`.
#' `visit_ids[[k]]` identifies those individuals in the exact row order of
#' `Y[[k]]`. Treatment is assigned once at baseline and shared across visits.
#'
#' @param X An `n` by `p` numeric matrix of baseline covariates (no missing
#'   values), with one row per individual.
#' @param T A binary numeric vector of length `n`, where 1 indicates the
#'   treatment of interest and 0 indicates the comparator.
#' @param Y A list of `m` numeric outcome matrices. The `k`th matrix has
#'   `n_k` rows (observed individuals at visit `k`) and `r` columns (outcome
#'   components). All matrices must have the same outcome columns and order.
#' @param visit_ids An optional list of `m` integer vectors. `visit_ids[[k]]`
#'   contains the row indices in `X` and `T` for individuals measured at
#'   visit `k`, in the same order as the rows of `Y[[k]]`. Individuals may
#'   occur at multiple visits, but indices within a visit must be unique.
#'   If `NULL`, all visits must contain all `n` individuals in baseline order.
#' @param eta A positive COMPACT ridge parameter. If `NULL` (the default),
#'   choose it by stratified subject-level cross-validation.
#' @param eta_grid Positive candidate ridge parameters used if `eta = NULL`.
#'   Default: `10^seq(-5, 1, length.out = 12)`.
#' @param nfolds Number of cross-validation folds used to select `eta`.
#'   Folds are defined at the individual level, so a subject never appears
#'   in both training and validation sets, even across different visits.
#'   Default: `5`.
#' @param seed Integer random seed for constructing cross-validation folds.
#'   Default: `1`.
#'
#' @details
#' At each visit, SCOUT constructs an initial COMPACT adjustment space for
#' treatment and all outcome components, then residualizes the observed
#' outcomes and treatment on that space. Its nonnegative least-squares
#' objective is the sum over visits of the residual sum of squares divided
#' by `n_k`, giving each visit equal influence regardless of sample size.
#' The NNLS solution is subsequently normalized to sum to one if nonzero;
#' `objective` is reported at the *unnormalized* NNLS optimum.
#' The final OLS uses the original COMPACT adjustment scores built for treatment
#' and all outcome components, with no phenotype-specific 2D refinement.
#'
#' The final score matrices contain predictions for **all** `n` baseline
#' individuals. ATEs and observed CATEs at visit `k` are evaluated only on
#' `visit_ids[[k]]`. The `CATE_all` output additionally contains score-based
#' predictions for all baseline individuals, including those not measured
#' at that visit.
#'
#' @return A list with components:
#' \describe{
#'   \item{alpha}{Length-`r` learned nonnegative outcome weights, normalized
#'     to sum to one when at least one weight is positive.}
#'   \item{active}{Indices of outcome components with weights above `1e-8`.}
#'   \item{objective}{Visit-normalized NNLS objective before rescaling weights.}
#'   \item{marginals}{Adjusted treatment-outcome associations standardized
#'     using pooled visit-normalized residual variation.}
#'   \item{marginals_by_visit}{An `r` by `m` matrix of analogous associations
#'     computed separately at each visit.}
#'   \item{eta}{Selected or user-supplied ridge parameter.}
#'   \item{eta_cv}{Cross-validation summary (or `NULL` when `eta` is supplied).}
#'   \item{eta_cv_folds}{Fold-by-parameter cross-validation scores, or `NULL`.}
#'   \item{learning_scores}{List of initial score matrices, each with `n`
#'     rows and at most `r + 1` columns.}
#'   \item{final}{List of visit-specific results. Each element contains
#'     `B` (initial COMPACT score directions, up to `r + 1`),
#'     `S` (`n`-row score matrix),
#'     `S_visit` (scores for observed individuals), `ids` (their indices),
#'     `ATE` (mean estimated effect among observed individuals),
#'     `CATE` (length-`n_k` vector of observed-subject effects),
#'     `CATE_all` (length-`n` predictions for all baseline individuals),
#'     and `coef` (OLS coefficients).}
#'   \item{visit_ids}{Validated index vectors defining each visit's sample.}
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
#' fit <- SCOUT_no_2D(X, T, Y, visit_ids = ids, eta = 0.1)
#' fit$alpha
#' fit$final[[2]]$ATE
#'
#' @seealso [SCOUT()]
#' @md
#' @export
SCOUT_no_2D <- function(X, T, Y, eta = NULL,
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
      nfolds != as.integer(nfolds) || nfolds < 2L)
    stop("nfolds must be an integer greater than or equal to two.")

  # Residuals and R^2 for the COMPACT cross-validation criterion.
  res <- function(A, Z) stats::lm.fit(cbind(1, Z), A)$residuals
  R2 <- function(y, Z) {
    Z <- as.matrix(Z)
    d <- sum((y - mean(y))^2)
    if (!ncol(Z) || d <= .Machine$double.eps) return(0)
    max(0, 1 - sum(res(y, Z)^2) / d)
  }

  # Precompute visit-specific COMPACT quantities. Each visit uses only
  # its observed individuals (and its own n_k and centering means).
  prep <- function(X, T, Y, visit_ids) {
    visits <- lapply(seq_along(Y), function(k) {
      ids <- visit_ids[[k]]
      Xk <- X[ids, , drop = FALSE]
      Tk <- T[ids]
      nk <- length(ids)
      xm <- colMeans(Xk)
      Xc <- sweep(Xk, 2, xm, "-")
      Tc <- Tk - mean(Tk)
      tt <- sum(Tc^2)
      if (tt <= 0)
        stop("Treatment does not vary within training visit ", k, ".")
      XX <- crossprod(Xc)
      xt <- crossprod(Xc, Tc)
      V <- (crossprod(Xc, Y[[k]]) -
            xt %*% (crossprod(Tc, Y[[k]]) / tt)) / nk
      G <- XX / nk - tcrossprod(drop(xt)) / (2 * nk * tt)
      E <- eigen((G + t(G)) / 2, symmetric = TRUE)
      E$values <- pmax(E$values, 0)
      list(X = Xc, T = Tk, Y = Y[[k]], ids = ids, n = nk,
           p = ncol(X), xm = xm, XX = XX, vT = xt / nk, V = V, E = E)
    })
    list(X = X, T = T, visits = visits, n = nrow(X), m = length(visits))
  }

  # Apply (G_k + eta * I)^(-1) using its eigendecomposition.
  solveG <- function(Pk, A, e)
    Pk$E$vectors %*%
      sweep(crossprod(Pk$E$vectors, A), 1, Pk$E$values + e, "/")

  # Whiten score directions in the observed visit sample. Different rotations
  # and rescalings preserve the same adjustment space.
  whiten <- function(Pk, W) {
    C <- crossprod(W, Pk$XX %*% W) / Pk$n
    E <- eigen((C + t(C)) / 2, symmetric = TRUE)
    keep <- E$values > 1e-10 * max(E$values)
    if (!any(keep)) return(matrix(0, Pk$p, 0L))
    W %*% sweep(E$vectors[, keep, drop = FALSE],
                2, sqrt(E$values[keep]), "/")
  }

  # Residualize on an intercept and whitened scores. Whitening makes
  # crossprod(S)/n_k the identity in each observed visit sample.
  adj <- function(A, S) {
    A <- as.matrix(A)
    A <- sweep(A, 2, colMeans(A), "-")
    if (ncol(S)) A <- A - S %*% (crossprod(S, A) / nrow(S))
    if (ncol(A) == 1L) drop(A) else A
  }

  fit <- function(P, e, effects = TRUE) {
    K <- P$m
    pk <- P$visits
    basis <- function(Pk, V) whiten(Pk, solveG(Pk, cbind(Pk$vT, V), e))

    # Stage 1: COMPACT scores for the treatment direction and all r outcomes.
    B0 <- lapply(pk, function(Pk) basis(Pk, Pk$V))
    S0_obs <- Map(function(Pk, B) Pk$X %*% B, pk, B0)
    RY <- Map(function(Pk, S) as.matrix(adj(Pk$Y, S)), pk, S0_obs)
    RT <- Map(function(Pk, S) as.numeric(adj(Pk$T, S)), pk, S0_obs)

    # Stage 2: learn one outcome weight vector shared across all visits.
    # Visit k contributes ||RT_k - RY_k %*% alpha||^2 / n_k.
    nk <- vapply(pk, function(Pk) Pk$n, 0L)
    H <- Reduce(`+`, Map(function(R, n) crossprod(R) / n, RY, nk))
    cvec <- Reduce(`+`, Map(function(R, t, n) drop(crossprod(R, t)) / n,
                            RY, RT, nk))
    tt <- sum(vapply(seq_len(K), function(k) sum(RT[[k]]^2) / nk[k], 0.0))
    if (tt <= .Machine$double.eps)
      stop("No residual treatment variation remains after initial adjustment.")
    fg <- function(a) tt - 2 * sum(a * cvec) + sum(a * drop(H %*% a))
    gr <- function(a) 2 * (drop(H %*% a) - cvec)
    o <- stats::optim(rep(0, r), fg, gr, method = "L-BFGS-B",
                      lower = rep(0, r),
                      control = list(maxit = 500, pgtol = 1e-8))
    a <- pmax(o$par, 0)
    obj <- fg(a)  # Objective evaluated before post-NNLS normalization.
    if (sum(a) > 0) a <- a / sum(a)

    # Optional descriptive adjusted marginal associations for the r outcomes.
    if (effects) {
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
    } else {
      marginals <- marginals_by_visit <- NULL
    }

    # Stage 3 (ablation): retain the initial COMPACT scores for final OLS.
    # Unlike SCOUT, do not refine to two phenotype-specific scores. The score
    # directions are learned from the observed visit sample and applied to all
    # baseline individuals; fit OLS only on the n_k observed outcomes.
    final <- lapply(seq_len(K), function(k) {
      Pk <- pk[[k]]
      B <- B0[[k]]
      if (!effects) return(list(B = B))

      S_all <- sweep(P$X, 2, Pk$xm, "-") %*% B
      S_obs <- S_all[Pk$ids, , drop = FALSE]
      y <- drop(Pk$Y %*% a)
      q <- ncol(S_obs)
      Z <- cbind(1, Pk$T, S_obs, S_obs * Pk$T)
      b <- stats::lm.fit(Z, y)$coefficients
      if (any(!is.finite(b)))
        stop("Rank-deficient final OLS model at visit ", k,
             "; check visit sample size and score variation.")
      cate <- if (q) drop(b[2] + S_obs %*% tail(b, q)) else rep(b[2], Pk$n)
      cate_all <- if (q) drop(b[2] + S_all %*% tail(b, q)) else rep(b[2], P$n)
      list(B = B, S = S_all, S_visit = S_obs, ids = Pk$ids,
           ATE = mean(cate), CATE = cate, CATE_all = cate_all, coef = b)
    })

    # Score matrices are returned for all baseline individuals.
    learning_scores <- if (effects) Map(function(Pk, B) {
      sweep(P$X, 2, Pk$xm, "-") %*% B
    }, pk, B0) else NULL

    list(alpha = a, objective = obj,
         learning_scores = learning_scores,
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

    cv <- vapply(eta_grid, function(e) vapply(folds, function(D) {
      M <- fit(D$P, e, effects = FALSE)
      mean(vapply(seq_len(m), function(k) {
        Pk <- D$P$visits[[k]]
        Xval <- sweep(D$Xval[[k]], 2, Pk$xm, "-")
        S <- Xval %*% M$final[[k]]$B
        if (!ncol(S)) return(0)
        y <- drop(D$Yval[[k]] %*% M$alpha)
        t <- D$Tval[[k]]
        sqrt(R2(t, S) * R2(res(y, t), res(S, t)))
      }, numeric(1)))
    }, numeric(1)), numeric(nfolds))

    mu <- colMeans(cv)
    eta <- eta_grid[which.max(mu)]
    eta_cv <- data.frame(eta = eta_grid, mean_score = mu,
                         se_score = apply(cv, 2, sd) / sqrt(nfolds))
  } else {
    cv <- eta_cv <- NULL
  }

  # Fit on the complete available data at each visit using the selected eta.
  M <- fit(prep(X, T, Y, visit_ids), eta, effects = TRUE)
  list(alpha = M$alpha,
       active = which(M$alpha > 1e-8),
       objective = M$objective,
       marginals = M$marginals,
       marginals_by_visit = M$marginals_by_visit,
       eta = eta,
       eta_cv = eta_cv,
       eta_cv_folds = cv,
       learning_scores = M$learning_scores,
       final = M$final,
       visit_ids = visit_ids)
}
