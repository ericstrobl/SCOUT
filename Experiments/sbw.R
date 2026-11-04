# Minimal SBW-SCORE comparator
# Learns alpha with weighted longitudinal SCORE, returns native SBW ATEs,
# then estimates tau(X) for the fixed phenotype with weighted ridge T x X.

.simplex <- function(x) {
  u <- sort(x, decreasing=TRUE); z <- cumsum(u)-1
  j <- max(which(u-z/seq_along(u)>0)); pmax(x-z[j]/j, 0)
}

.wcor <- function(x, y, w) {
  w <- w/sum(w); x <- x-sum(w*x); y <- y-sum(w*y)
  d <- sqrt(sum(w*x^2)*sum(w*y^2))
  if (!is.finite(d) || d < 1e-14) 0 else sum(w*x*y)/d
}

.score_fit <- function(T, Y, w, starts=10, maxit=1000, tol=1e-8, seed=1) {
  obj <- function(a) sum(vapply(Y, \(M) .wcor(T, drop(M%*%a), w), 0.0))
  grad <- function(a) {
    ww <- w/sum(w); t <- T-sum(ww*T); nt <- sqrt(sum(ww*t^2)); g <- 0*a
    if (nt < 1e-14) return(g)
    for (M in Y) {
      Mc <- sweep(M, 2, colSums(ww*M), "-"); z <- drop(Mc%*%a); nz2 <- sum(ww*z^2)
      if (nz2 < 1e-24) next
      nz <- sqrt(nz2); rho <- sum(ww*t*z)/(nt*nz)
      g <- g + drop(crossprod(Mc,ww*t))/(nt*nz) - rho*drop(crossprod(Mc,ww*z))/nz2
    }
    g
  }
  r <- ncol(Y[[1]]); set.seed(seed)
  item <- vapply(seq_len(r), \(j) {a<-numeric(r); a[j]<-1; obj(a)}, 0.0)
  A <- list(rep(1/r,r), replace(numeric(r), which.max(item), 1))
  if (starts > 2) for (j in 3:starts) {a<-rexp(r); A[[j]]<-a/sum(a)}
  one <- function(a) {
    a <- .simplex(a); q <- obj(a)
    for (it in seq_len(maxit)) {
      g <- grad(a); step <- 1; ok <- FALSE
      for (ls in 1:50) {
        b <- .simplex(a+step*g); qb <- obj(b)
        if (is.finite(qb) && qb >= q + 1e-4*sum(g*(b-a)) - 1e-12) {ok<-TRUE; break}
        step <- step/2
      }
      if (!ok || sum(abs(b-a)) < tol || abs(qb-q) < tol) {if(ok){a<-b;q<-qb}; break}
      a <- b; q <- qb
    }
    list(alpha=a, objective=q)
  }
  z <- lapply(A[seq_len(min(starts,length(A)))], one)
  z[[which.max(vapply(z, `[[`, 0.0, "objective"))]]
}

.ridge_cate <- function(X,T,y,w,nfolds=5,lambda="lambda.min",seed=1) {
  p <- ncol(X); K <- min(nfolds,min(table(T))); set.seed(seed)
  foldid <- integer(length(T)); for(g in 0:1){ii<-which(T==g); foldid[ii]<-sample(rep(seq_len(K),length.out=length(ii)))}
  Z <- cbind(T,X,X*T)
  cv <- glmnet::cv.glmnet(Z,y,alpha=0,weights=w,foldid=foldid,type.measure="mse",
                          penalty.factor=c(0,rep(1,2*p)),standardize=TRUE)
  drop(predict(cv,cbind(1,X,X),s=lambda) - predict(cv,cbind(0,X,0*X),s=lambda))
}

sbw_score <- function(X,T,Y_list,balance_tol=.05,target_tol=.05,solver="osqp",
                      score_starts=10,score_maxit=1000,score_tol=1e-8,
                      cate_nfolds=5,cate_lambda="lambda.min",B=0,seed=1,
                      verbose=TRUE,...) {
  t0 <- proc.time()[["elapsed"]]
  if (B > 0) stop("Minimal sbw_score() is point-estimation only; use B=0.")
  X <- as.matrix(X); T <- as.numeric(T); Y_list <- lapply(Y_list,as.matrix)
  d <- as.data.frame(X); names(d) <- make.names(if(is.null(colnames(X))) paste0("X",seq_len(ncol(X))) else colnames(X),unique=TRUE); d$.T <- T
  wf <- optweight::optweight(.T~.,data=d,estimand="ATE",tols=balance_tol,target.tols=target_tol,
                             norm="l2",std.cont=TRUE,min.w=1e-8,solver=solver,verbose=FALSE)
  w <- as.numeric(wf$weights)
  sf <- .score_fit(T,Y_list,w,score_starts,score_maxit,score_tol,seed)
  a <- sf$alpha; if(!is.null(colnames(Y_list[[1]]))) names(a) <- colnames(Y_list[[1]])
  ys <- lapply(Y_list, \(Y) drop(Y%*%a))
  wmean <- function(y,g) sum(w[T==g]*y[T==g])/sum(w[T==g])
  tau <- vapply(ys, \(y) wmean(y,1)-wmean(y,0), 0.0)
  CATE <- lapply(seq_along(ys), \(k) .ridge_cate(X,T,ys[[k]],w,cate_nfolds,cate_lambda,seed+200000+k))
  list(alpha=a,tau=tau,CATE=CATE,weights=w,stage1_objective=sf$objective,
       effective_sample_size=c(control=sum(w[T==0])^2/sum(w[T==0]^2),treated=sum(w[T==1])^2/sum(w[T==1]^2)),
       runtime_sec=unname(proc.time()[["elapsed"]]-t0))
}
