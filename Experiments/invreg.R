# Inverse-regression comparator
# Learns alpha at the last LOCF visit
# Estimates tau(X) for the fixed phenotype with IPW-weighted ridge T x X.

.ridge_cate <- function(X,T,y,w,nfolds=5,lambda="lambda.min",seed=1) {
  p <- ncol(X); K <- min(nfolds,min(table(T))); set.seed(seed)
  foldid <- integer(length(T)); for(g in 0:1){ii<-which(T==g); foldid[ii]<-sample(rep(seq_len(K),length.out=length(ii)))}
  Z <- cbind(T,X,X*T)
  cv <- glmnet::cv.glmnet(Z,y,alpha=0,weights=w,foldid=foldid,type.measure="mse",
                          penalty.factor=c(0,rep(1,2*p)),standardize=TRUE)
  drop(predict(cv,cbind(1,X,X),s=lambda) - predict(cv,cbind(0,X,0*X),s=lambda))
}

invreg_nn <- function(X,T,Y_list,nfolds=5,ps_clip=1e-4,lambda="lambda.min",
                      cate_nfolds=5,cate_lambda="lambda.min",B=0,seed=1,
                      verbose=TRUE,...) {
  t0 <- proc.time()[["elapsed"]]
  if (B > 0) stop("Minimal invreg_nn() is point-estimation only; use B=0.")
  X <- as.matrix(X); T <- as.numeric(T); Y_list <- lapply(Y_list,as.matrix); n <- nrow(X)
  
  # Item-wise LOCF, then learn the inverse-regression phenotype at the last visit.
  Y <- Y_list[[1]]
  if(length(Y_list)>1) for(k in 2:length(Y_list)){Z<-Y_list[[k]]; miss<-is.na(Z); Z[miss]<-Y[miss]; Y<-Z}
  keep <- complete.cases(X,T,Y); X0 <- X[keep,,drop=FALSE]; T0 <- T[keep]; Y <- Y[keep,,drop=FALSE]
  
  # Ridge propensity and ATE IPW.
  K <- min(nfolds,min(table(T0))); set.seed(seed)
  foldid <- integer(length(T0)); for(g in 0:1){ii<-which(T0==g); foldid[ii]<-sample(rep(seq_len(K),length.out=length(ii)))}
  ps <- glmnet::cv.glmnet(X0,T0,family="binomial",alpha=0,foldid=foldid,type.measure="deviance",standardize=TRUE)
  e <- drop(predict(ps,X0,s=lambda,type="response")); e <- pmin(pmax(e,ps_clip),1-ps_clip)
  w <- T0/e + (1-T0)/(1-e)
  
  # Marginal item effects and nonnegative inverse regression T ~ Y.
  D <- cbind(1,T0)
  tau <- solve(crossprod(D,D*w),crossprod(D,Y*w))[2,]
  ww <- w/sum(w); Yc <- sweep(Y,2,colSums(ww*Y),"-"); Tc <- T0-sum(ww*T0); sw <- sqrt(w)
  beta <- as.numeric(coef(nnls::nnls(Yc*sw,sw*Tc))); names(beta) <- colnames(Y)
  sb <- sum(beta); ok <- is.finite(sb) && sb > 1e-12
  a <- if(ok) beta/sb else rep(NA_real_,length(beta)); names(a) <- names(beta)
  
  # CATE for the fixed phenotype at every visit, conditional on X.
  CATE <- if(ok) lapply(seq_along(Y_list), \(k) {
    y <- drop(Y_list[[k]][keep,,drop=FALSE]%*%a)
    .ridge_cate(X0,T0,y,w,cate_nfolds,cate_lambda,seed+200000+k)
  }) else lapply(Y_list, \(Y) rep(NA_real_,sum(keep)))
  
  list(alpha=a,alpha_defined=ok,beta=beta,ATE_alpha=if(ok) sum(a*tau) else NA_real_,
       ATEc=sum(beta*tau),marginal_ATE=tau,propensity=e,ipw=w,CATE=CATE,
       used_rows=which(keep),runtime_sec=unname(proc.time()[["elapsed"]]-t0))
}
