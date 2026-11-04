# Compact observational CV-TMLE comparator used by SCOUT simulations.
# Requires glmnet + quadprog.

.folds <- function(A, V, seed) {
  set.seed(seed); g <- split(seq_along(A), A); V <- min(V, min(lengths(g)))
  if (V < 2) stop("Too few observations per treatment group.")
  f <- integer(length(A)); for (ii in g) f[ii] <- sample(rep(seq_len(V), length.out=length(ii))); f
}

.pd <- function(S, eps=1e-8) {
  S <- (S+t(S))/2; e <- eigen(S, symmetric=TRUE)
  e$values <- pmax(e$values, max(eps*max(abs(e$values),1), eps))
  e$vectors %*% (e$values * t(e$vectors))
}

.opt_alpha <- function(mu, S) {
  K <- length(mu); sd0 <- sqrt(diag(S))
  if (all(mu <= 0)) { a <- numeric(K); a[which.max(mu/sd0)] <- 1 }
  else {
    q <- try(quadprog::solve.QP(S, mu, diag(K), rep(0,K))$solution, silent=TRUE)
    if (inherits(q,"try-error") || sum(pmax(q,0)) <= 1e-12) {
      a <- numeric(K); a[which.max(mu/sd0)] <- 1
    } else { a <- pmax(q,0); a <- a/sum(a) }
  }
  list(alpha=a, stat=sum(a*mu)/sqrt(drop(crossprod(a,S%*%a))))
}

.mc_p <- function(Tobs, S, B, seed) {
  set.seed(seed); K <- ncol(S); e <- eigen(S, symmetric=TRUE)
  Z <- matrix(rnorm(B*K),B,K) %*% t(e$vectors %*% diag(sqrt(pmax(e$values,0)),K))
  s <- apply(Z,1,function(z) .opt_alpha(z,S)$stat)
  mean(s >= Tobs)
}

.fit_nuisance <- function(X, A, Y, nfolds, seed, lambda_rule) {
  X <- as.matrix(X); Y <- as.matrix(Y); Z <- cbind(A, X, X*A)
  K <- min(nfolds, nrow(X)); if (K < 3) stop("Need >=3 observations for nuisance CV.")
  q <- lapply(seq_len(ncol(Y)), function(j) {
    set.seed(seed+j); foldid <- sample(rep(seq_len(K), length.out=nrow(X)))
    glmnet::cv.glmnet(Z,Y[,j],family="gaussian",alpha=0,foldid=foldid,
                      standardize=TRUE,penalty.factor=c(0,rep(1,2*ncol(X))))
  })
  gf <- .folds(A,nfolds,seed+50000)
  g <- glmnet::cv.glmnet(X,A,family="binomial",alpha=0,foldid=gf,standardize=TRUE)
  list(q=q,g=g,s=lambda_rule)
}

.pred_nuisance <- function(f, X, A, clip=.01) {
  X <- as.matrix(X); n <- nrow(X)
  pred <- function(z,M) drop(predict(z,newx=M,s=f$s,type="response"))
  g <- pmin(pmax(pred(f$g,X),clip),1-clip)
  Za <- cbind(A,X,X*A); Z1 <- cbind(1,X,X); Z0 <- cbind(0,X,0*X)
  getq <- function(M) sapply(f$q, pred, M=M)
  qo <- getq(Za); q1 <- getq(Z1); q0 <- getq(Z0)
  if (is.null(dim(qo))) { qo <- matrix(qo,ncol=1); q1 <- matrix(q1,ncol=1); q0 <- matrix(q0,ncol=1) }
  list(g=g,qo=qo,q1=q1,q0=q0)
}

.tmle <- function(y,A,qo,q1,q0,g) {
  H <- A/g-(1-A)/(1-g); den <- sum(H^2)
  eps <- if (den>0) sum(H*(y-qo))/den else 0
  qo <- qo+eps*H; q1s <- q1+eps/g; q0s <- q0-eps/(1-g)
  psi <- mean(q1s-q0s); ic <- H*(y-qo)+q1s-q0s-psi
  list(psi=psi,ic=ic)
}

.endpoint_tmle <- function(Y,A,p) {
  z <- lapply(seq_len(ncol(Y)), function(j) .tmle(Y[,j],A,p$qo[,j],p$q1[,j],p$q0[,j],p$g))
  list(psi=vapply(z,`[[`,0.0,"psi"), IC=do.call(cbind,lapply(z,`[[`,"ic")))
}

.prepare_Y <- function(Y_list, mode, higher_better) {
  if (match.arg(mode,c("last","stack")) == "last") Y <- as.matrix(Y_list[[length(Y_list)]])
  else Y <- do.call(cbind,lapply(Y_list,as.matrix))
  if (!higher_better) Y <- -Y
  Y
}

cv_tmle_global <- function(X,T,Y_list,treated=1,higher_better=TRUE,
                           outcome_mode=c("last","stack"),V=10,C=.25,B_null=5000,
                           alpha_ref=NULL,nuisance_folds=5,lambda_rule="lambda.min",
                           g_clip=.01,cov_epsilon=1e-8,seed=1,verbose=TRUE) {
  t0 <- proc.time()[["elapsed"]]; X <- as.matrix(X); A <- as.numeric(T==treated)
  Y <- .prepare_Y(Y_list,outcome_mode,higher_better); n <- nrow(X); K <- ncol(Y)
  if (is.null(alpha_ref)) alpha_ref <- rep(1/K,K)
  folds <- .folds(A,V,seed); V <- max(folds)
  aa <- ast <- matrix(NA_real_,V,K); psi <- numeric(V); ic <- numeric(n)
  for (v in seq_len(V)) {
    if (verbose) message("CV-TMLE fold ",v,"/",V)
    tr <- folds!=v; va <- !tr
    nf <- .fit_nuisance(X[tr,,drop=FALSE],A[tr],Y[tr,,drop=FALSE],nuisance_folds,seed+100*v,lambda_rule)
    pt <- .pred_nuisance(nf,X[tr,,drop=FALSE],A[tr],g_clip)
    et <- .endpoint_tmle(Y[tr,,drop=FALSE],A[tr],pt); S <- .pd(cov(et$IC),cov_epsilon)
    op <- .opt_alpha(et$psi,S); aa[v,] <- op$alpha
    p <- .mc_p(sqrt(sum(tr))*op$stat,S,B_null,seed+10000+v)
    sh <- min(1,C*log(n)*p); ast[v,] <- (1-sh)*aa[v,]+sh*alpha_ref
    pv <- .pred_nuisance(nf,X[va,,drop=FALSE],A[va],g_clip); a <- ast[v,]
    z <- .tmle(drop(Y[va,,drop=FALSE]%*%a),A[va],drop(pv$qo%*%a),drop(pv$q1%*%a),drop(pv$q0%*%a),pv$g)
    psi[v] <- z$psi; ic[va] <- z$ic
  }
  est <- mean(psi); se <- sd(ic)/sqrt(n); alpha_mean <- colMeans(ast)
  list(estimate=est,se=se,alpha_mean=alpha_mean,alpha=alpha_mean,
       alpha_stabilized=ast,alpha_adaptive=aa,folds=folds,fold_estimates=psi,
       runtime_sec=proc.time()[["elapsed"]]-t0)
}

tmle_fixed_alpha <- function(X,T,Y,alpha,treated=1,seed=1,V=10,nuisance_folds=5,
                             lambda_rule="lambda.min",g_clip=.01,return_cate=FALSE) {
  X <- as.matrix(X); A <- as.numeric(T==treated); y <- drop(as.matrix(Y)%*%alpha)
  folds <- .folds(A,V,seed); psi <- numeric(max(folds)); cate <- numeric(nrow(X))
  for (v in seq_len(max(folds))) {
    tr <- folds!=v; va <- !tr
    nf <- .fit_nuisance(X[tr,,drop=FALSE],A[tr],matrix(y[tr],ncol=1),nuisance_folds,seed+100*v,lambda_rule)
    p <- .pred_nuisance(nf,X[va,,drop=FALSE],A[va],g_clip)
    cate[va] <- drop(p$q1-p$q0)  # cross-fitted plug-in CATE; TMLE fluctuation targets ATE
    psi[v] <- .tmle(y[va],A[va],drop(p$qo),drop(p$q1),drop(p$q0),p$g)$psi
  }
  if (return_cate) list(ATE=mean(psi),CATE=cate,folds=folds) else mean(psi)
}