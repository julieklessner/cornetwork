#' SparCC
#'
#' reimplementation of SparCC in R
#'
#' @param data count data
#' @param iter number of iterations to compute median correlation
#' @param inner_iter number of iterations for the inner loop
#' @export
sparcc <- function(data, iter=20, inner_iter=10, th=.1) {
  ## 
  #  without all the 'frills'
  sparccs <- 
    lapply(1:iter, function(i) 
      sparccinner(t(apply(data, 1, norm_diric)), iter=inner_iter, th=th))
  # collect
  cors <- array(unlist(lapply(sparccs, function(x) x$Cor)),
                c(ncol(data),ncol(data),iter))
  corMed <- apply(cors, 1:2, median)
  covs <- array(unlist(lapply(sparccs, function(x) x$Cov)),
                c(ncol(data),ncol(data),iter))
  covMed <- apply(cors, 1:2, median)
  
  covMed <- cor2cov(corMed, sqrt(diag(covMed)))
  list(Cov=covMed, Cor=corMed)
}


#' @keywords internal
sparccinner <- function(data.f, T=NULL, iter=10, th=0.1) {
  if (is.null(T))   T  <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor
  
  ## do iterations here
  excluded <- NULL
  for (i in 1:iter) {
    res.excl <- exclude_pairs(Cor, M, th, excluded)
    M <- res.excl$M
    excluded <- res.excl$excluded
    if (res.excl$break_flag) break
    res.bv <- basis_var(T, M=M, excluded=excluded)
    Vbase  <- res.bv$Vbase
    M      <- res.bv$M
    K <- M
    diag(K) <- 1
    cbase  <- C_from_V(T, Vbase)
    Cov    <- cbase$Cov
    Cor    <- cbase$Cor
  }
  list(Cov=Cov, Cor=Cor, i=i, M=M, excluded=excluded)
}

#' @keywords internal
exclude_pairs <- function(Cor, M, th=0.1, excluded=NULL) {
  # exclude pairs with high correlations
  break_flag <- FALSE
  C_temp <- abs(Cor - diag(diag(Cor)) )  # abs value / remove diagonal
  if (!is.null(excluded)) C_temp[excluded] <- 0 # set previously excluded correlations to 0
  exclude <- which(abs(C_temp - max(C_temp)) < .Machine$double.eps*100)[1:2]
  if (max(C_temp) > th)  {
    i <- na.exclude(arrayInd(exclude, c(nrow(M), ncol(M)))[,1])
    M[i,i] <- M[i,i] - 1
    excluded_new <- c(excluded, exclude)
  } else {
    excluded_new <- excluded
    break_flag   <- TRUE
  }
  list(M=M, excluded=excluded_new, break_flag=break_flag)
}

#' @keywords internal
basis_cov <- function(data.f) {
  # data.f -> relative abundance data
  # OTUs in columns, samples in rows (yes, I know this is transpose of normal)
  # first compute aitchison variation
  T <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor
  list(Cov=Cov, M=M)
}

#' @keywords internal
basis_var <- function(T, CovMat = matrix(0, nrow(T), ncol(T)), 
                      M = matrix(1, nrow(T), ncol(T)) + (diag(ncol(T))*(ncol(T)-2)), 
                      excluded = NULL, Vmin=1e-4) {
  
  if (!is.null(excluded)) {
    T[excluded] <- 0
    #   CovMat[excluded] <- 0
  }
  Ti     <- matrix(rowSums(T))
  CovVec <- matrix(rowSums(CovMat - diag(diag(CovMat)))) # row sum of off diagonals
  M.I <- tryCatch(solve(M), error=function(e) MASS::ginv(M))
  Vbase <- M.I %*% (Ti + 2*CovVec)
  Vbase[Vbase < Vmin] <- Vmin
  list(Vbase=Vbase, M=M)
}

#' @keywords internal
C_from_V <- function(T, Vbase) {
  J      <- matrix(1, nrow(T), ncol(T))
  Vdiag  <- diag(c(Vbase))
  CovMat <- .5*((J %*% Vdiag) + (Vdiag %*% J) - T)
  CovMat <- (CovMat + t(CovMat))/2  # enforce symmetry
  # check that correlations are within -1,1
  CorMat <- cov2cor(CovMat)
  CorMat[abs(CorMat) > 1] <- sign(CorMat[abs(CorMat) > 1])
  CovMat <- cor2cov(CorMat, sqrt(as.vector(Vbase)))
  list(Cov=CovMat, Cor=CorMat)
}

#' @keywords internal
av <- function(data) {
  cov.clr <- cov(clr(data))
  J <- matrix(1, ncol(data), ncol(data))
  (J %*% diag(diag(cov.clr))) + (diag(diag(cov.clr)) %*% J) - (2*cov.clr)
}


#' importFrom VGAM rdiric
#' @export
norm_diric   <- function(x, rep=1) {
  dmat <- VGAM::rdiric(rep, x+1)
  norm_to_total(colMeans(dmat))
}

#####################################################################
# Different normalization schemes for microbiome counts (real or fake)
#
# @author Zachary Kurtz
# @date 10/10/2013
#####################################################################

#' add pseudocount before normalizing a vector
#' @export
norm_pseudo  <- function(x) norm_to_total(x+1)


#' Total sum normalization of a (presumably positive) data vector
#' @export
norm_to_total <- function(x) x/sum(x)


#' compute the shannon entropy from a vector (normalized internally)
#'
#' Shannon entropy is:
#'     sum [ x_i log(x_i) ]
#'      
#' @param x data vector
#' @return shannon entropy in base e
#' @export
shannon <- function(x) {
  x.f <- norm_to_total(x)
  -sum(x.f*log(x.f), na.rm=TRUE)
}

#' N_effective: Compute the exponential of the shannon entropy. linearizes shannon entropy, for a better diveristy metric (effective number of species)
#'
#' @param x data vector
#' @return N_eff in base e
#' @export
neff <- function(x) exp(shannon(x))


#' Centered log-ratio functions
#' @export
clr <- function(x, ...) {
  UseMethod('clr', x)
}

#' @method clr default
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
  apply(x.f, mar, clr, ...)
}

#' @method clr data.frame
#' @export
clr.data.frame <- function(x.f, mar=2, ...) {
  clr(as.matrix(x.f), mar, ...)
}

#' Additive log-ratio functions
#' @export
alr <- function(x, ...) {
  UseMethod("alr", x)
}

#' @method alr default
#' @export
alr.default <- function(x.f, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps) {
  zero <- (x.f >= tol)
  LOG <- log(ifelse(zero, x.f, 1), base)
  x.alr <- ifelse(zero, LOG - LOG[divcomp], 0.0)
  if (removeDivComp) x.alr[-divcomp]
  else x.alr
}


#' @method alr matrix
#' @export
alr.matrix <- function(x.f, mar=2, divcomp=1, base=exp(1), removeDivComp=TRUE,
                       tol=.Machine$double.eps) {
  if (mar == 1) x.f <- t(x.f)
  zero <- (x.f >= tol)
  LOG <- log(ifelse(zero, x.f, 1), base)
  x.alr <- ifelse(zero, LOG - LOG[,divcomp], 0.0)
  if (removeDivComp) x.alr[,-divcomp]
  else x.alr
}

#' @method alr data.frame
#' @export
alr.data.frame <- function(x.f, mar=2, ...) {
  alr(as.matrix(x.f), mar, ...)
}

#' @export
triu <- function(x) x[upper.tri(x)]
#' @export
tril <- function(x) x[lower.tri(x)]

#' @export
triu2diag <- function(x, diagval=0) {
  e <- length(x)
  n <- .5 * (sqrt(8*e + 1)+1)
  mat <- matrix(0, n, n)
  mat[upper.tri(mat)] <- x
  mat <- mat + t(mat)
  diag(mat) <- diagval
  mat
}

#' Draw samples from a zero-inflated poisson distribution
#'
#' @param n the number of samples to draw
#' @param lambda The poisson rate parameter
#' @param pstr0 probability of drawing a zero
#' @return Poisson counts of length \eqn{n}
#' @importFrom stats qpois dpois runif
#' @export
rzipois <- function(n, lambda, pstr0 = 0) {
  ans <- rpois(n, lambda)
  ans <- ifelse(runif(n) < pstr0, 0, ans)
  prob0 <- exp(-lambda)
  deflat.limit <- -1/expm1(lambda)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * pstr0[ind0]
    ans[ind0] <- qpois(p = runif(sum(ind0), 
                                 min = dpois(0, lambda[ind0])), lambda[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN
  ans
}


#' @keywords internal
.zipois_getLam <- function(mu, S) {
  S <- max(sqrt(mu), S)
  (S^2/mu) + mu - 1
}

#' @keywords internal
.zipois_getP <- function(mu, S) {
  S <- max(sqrt(mu), S)
  (S^2 - mu) / (mu^2 - mu + S^2)
}

#' Generate multivariate, Zero-inflated poisson data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param lambdas supply rate parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzipois
#' @export
rmvzipois <- function(n, mu, Sigma, lambdas, ps, ...) {
  d   <- ncol(Sigma)
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  
  if (missing(lambdas) || missing(ps)) {
    if (missing(mu)) stop("Need to supply mu")
    if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
    lambdas <- unlist(lapply(1:length(SDs), function(i) .zipois_getLam(mu[i], SDs[i])))
    ps   <- unlist(lapply(1:length(SDs), function(i) .zipois_getP(mu[i], SDs[i])))
  }
  if (length(lambdas) != length(SDs)) stop("Sigma and mu/lambdas dimensions don't match")
  if (length(lambdas) == 1) stop("Need more than 1 variable")
  
  normd  <- rmvnorm(n, rep(0, d), Cor)
  unif   <- pnorm(normd)
  data <- t(apply(unif, 1, function(probs) {
    VGAM::qzipois(probs, lambdas, pstr0=ps, ...)
  }))
  
  
  data <- .fixInf(data)
  return(data)
}


#' Generate multivariate poisson data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param lambdas supply rate parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom stats qpois
#' @export
rmvpois <- function(n, mu, Sigma, ...) {
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  d   <- length(SDs)
  if (length(mu) != length(SDs)) stop("Sigma and mu/lambdas dimensions don't match")
  if (length(mu) == 1) stop("Need more than 1 variable")
  normd  <- rmvnorm(n, rep(0, d), Cor)
  unif   <- pnorm(normd)
  data <- t(qpois(t(unif), mu, ...))
  data <- .fixInf(data)
  return(data)
}

#' @keywords internal
.negbin_getK <- function(mu, S) {
  S <- max(sqrt(mu), S)
  mu^2/((S^2+1e-3)-mu)
}

#' Generate multivariate, Zero-inflated negative binomial data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param munbs Rate/mean parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @param ks shape parameter
#' @param ... other arguments to the negative binomial distribution
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom stats qnbinom
#' @export
rmvnegbin <- function(n, mu, Sigma, ks, ...) {
  # Generate an NxD matrix of Zero-inflated poisson data,
  # with counts approximately correlated according to Sigma
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(mu)) stop('mu is required')
  if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
  if (missing(ks)) {
    ks   <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i], SDs[i])))
  }
  d   <- length(mu)
  normd  <- rmvnorm(n, rep(0, d), Sigma=Cor)
  unif   <- pnorm(normd)
  data <- t(qnbinom(t(unif), mu=mu, size=ks, ...))
  data <- .fixInf(data)
  return(data)
}


#' @keywords internal
.zinegbin_getLam <- function(mu, S) {
  S   <- max(sqrt(mu)+1e-3, S)
  (mu + (mu^2 - mu + S^2) / mu) / 2
}

#' @keywords internal
.zinegbin_getP <- function(mu, lam) {
  1 - (mu / lam)
}

#' @keywords internal
.zinegbin_getK <- function(mu, S, lam) {
  S   <- max(sqrt(mu)+1e-3, S)
  (mu * lam) / (mu^2 - (mu * (lam + 1)) + S^2)
}


#' Generate multivariate, negative binomial data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param munbs Rate/mean parameter (instead of mu)
#' @param ks shape parameter
#' @param ... other arguments to the negative binomial distribution
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzinegbin
#' @export
rmvzinegbin <- function(n, mu, Sigma, munbs, ks, ps, ...) {
  # Generate an NxD matrix of Zero-inflated poisson data,
  # with counts approximately correlated according to Sigma
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(munbs) || missing(ps) || missing(ks)) {
    if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
    munbs <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getLam(mu[i], SDs[i])))
    ps   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getP(mu[i], munbs[i])))
    ks   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getK(mu[i], SDs[i], munbs[i])))
  }
  if (length(munbs) != length(SDs)) stop("Sigma and mu dimensions don't match")
  d   <- length(munbs)
  normd  <- rmvnorm(n, rep(0, d), Sigma=Cor)
  unif   <- pnorm(normd)
  data <- t(apply(unif, 1, function(probs) {
    VGAM::qzinegbin(probs, munb=munbs, size=ks, pstr0=ps, ...)
  }))
  data <- .fixInf(data)
  return(data)
}



#' @keywords internal
.fixInf <- function(data) {
  # hacky way of replacing infinite values with the col max + 1
  if (any(is.infinite(data))) {
    data <-  apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind<-which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm=TRUE)+1
      }
      x
    })
  }
  data
}


#' Draw samples from multivariate, correlated normal distribution
#' with counts correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param tol numerical tolerance for a zero eigenvalue (check for PD of Sigma)
#' @param empirical is Sigma the empirical correlation?
#' @return \eqn{Dxn} matrix with Gaussian data
#' @export
rmvnorm <- function(n=100, mu=rep(0,10), Sigma=diag(10), tol=1e-6, empirical=TRUE) {
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0, nv = length(mu))$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  return(t(X))
}

#' @export
multivnomial <- function(n, mu, Sigma) {
  N <- sum(mu)
  data <- rmultinom(n, N, mu/N)
  c <- as.matrix(chol(Sigma))
  t(data) %*% c
}

#' Convert a symmetric correlation matrix to a covariance matrix
#' given the standard deviation
#'
#' @param cor a symmetric correlation matrix
#' @param sds standard deviations of the resulting covariance.
#' @return Covariance matrix of sample dimension as cor
#' @export
cor2cov <- function(cor, sds) {
  if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
  cor * sds * rep(sds, each=nrow(cor))
}

