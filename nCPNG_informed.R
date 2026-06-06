## ============================================================================
##  Locally-balanced informed sampler for Bayesian Gaussian DAG learning
##  under the non-conjugate Normal-Gamma (nCPNG) prior.
##
##  Companion code for the JCGS resubmission
##  "Informed and Scalable Sampling for Bayesian Gaussian DAG Learning
##   under Non-conjugate Normal-Gamma Priors"
##  (S. Nazari, M. Arashi, A. Sadeghkhani)
##
##  This file implements:
##    (1) the closed-form node score of Theorem 1 (log-stable),
##        with an exact Bessel evaluation and a Bessel-free uniform
##        asymptotic approximation (Section 5);
##    (2) the random-walk Metropolis-Hastings sampler (baseline);
##    (3) the locally-balanced informed sampler of Algorithm 1,
##        with incremental neighbour-score caching;
##    (4) effective-sample-size and timing comparisons.
##
##  Repository: https://github.com/Samaneh-Nazari/nCPNG-DAG
## ============================================================================

suppressPackageStartupMessages({
  for (p in c("mvtnorm","huge","igraph","gRbase","GIGrvg","Bessel","coda"))
    if (!requireNamespace(p, quietly = TRUE))
      message(sprintf("Install package '%s' before running.", p))
})
library(mvtnorm); library(gRbase); library(GIGrvg); library(Bessel); library(coda)

## ----------------------------------------------------------------------------
## 1.  Closed-form node score (Theorem 1), log scale
## ----------------------------------------------------------------------------

## Exact log Bessel via Bessel::BesselK with exponential scaling for stability:
##   log K_nu(z) = log( BesselK(z, nu, expon.scaled=TRUE) ) - z
log_besselK_exact <- function(z, nu) {
  log(BesselK(z, nu, expon.scaled = TRUE)) - z
}

## Bessel-free uniform asymptotic approximation (DLMF 10.41), used when |nu|
## is large.  Returns log K_nu(z).  Includes the u_1 correction term so that
## the relative error is O(nu^{-2}) (Proposition 7 in the paper).
log_besselK_uniform <- function(z, nu) {
  nu <- abs(nu)
  zz <- z / nu                       # argument scaled so that z = nu * zz
  s  <- sqrt(1 + zz^2)
  t  <- 1 / s
  eta <- s + log(zz / (1 + s))
  ## leading term:  sqrt(pi/(2 nu)) * exp(-nu eta) / (1+z'^2)^{1/4}
  log_lead <- 0.5 * log(pi / (2 * nu)) - nu * eta - 0.25 * log(1 + zz^2)
  ## first correction u_1(t) = (3t - 5t^3)/24 ; factor (1 + u_1/nu)
  u1 <- (3 * t - 5 * t^3) / 24
  log_lead + log1p(u1 / nu)
}

## Switch: use uniform asymptotics when the order is large enough
log_besselK <- function(z, nu, approx = FALSE, thresh = 30) {
  if (approx && abs(nu) >= thresh) log_besselK_uniform(z, nu)
  else                              log_besselK_exact(z, nu)
}

pa_set <- function(node, DAG) which(DAG[, node] != 0)

## Node score:  log m(X_j | X_{pa(j)}, D)   (Theorem 1)
node_score <- function(node, DAG, tXX, n, a, U, approx = FALSE) {
  j <- node; pj <- pa_set(j, DAG); q <- ncol(tXX)
  a.star <- a + length(pj) - q + 1
  Upost  <- U + tXX
  if (length(pj) == 0) {
    U_jj  <- U[j, j]
    kappa <- 0.5 * tXX[j, j]
    lambda<- 0.5 * U_jj
    nu    <- 0.5 * (a.star - n)
    zarg  <- 2 * sqrt(kappa * lambda)
    val <- -n/2*log(2*pi) + a.star/2*log(U_jj/2) - lgamma(a.star/2) +
           log(2) + log_besselK(zarg, nu, approx) +
           (nu/2) * (log(kappa) - log(lambda))
  } else {
    U_pp   <- U[pj, pj, drop = FALSE]
    Up_pp  <- Upost[pj, pj, drop = FALSE]
    iU_pp  <- solve(U_pp)
    iUp_pp <- solve(Up_pp)
    U_paj  <- U[pj, j]; Up_paj <- Upost[pj, j]
    U_jj   <- U[j, j] - t(U_paj) %*% iU_pp %*% U_paj
    kappa  <- as.numeric(0.5 * (-t(Up_paj) %*% iUp_pp %*% Up_paj +
                                 t(U_paj)  %*% iU_pp  %*% U_paj + tXX[j, j]))
    lambda <- as.numeric(0.5 * U_jj)
    nu     <- 0.5 * (a.star - n)
    zarg   <- 2 * sqrt(max(kappa, 1e-12) * max(lambda, 1e-12))
    logdet <- determinant(U_pp, logarithm = TRUE)$modulus -
              determinant(Up_pp, logarithm = TRUE)$modulus
    val <- -n/2*log(2*pi) + 0.5*logdet + a.star/2*log(U_jj/2) -
           lgamma(a.star/2) + log(2) + log_besselK(zarg, nu, approx) +
           (nu/2) * (log(kappa) - log(lambda))
  }
  if (!is.finite(as.numeric(val))) return(-1e10)
  Re(as.numeric(val))
}

## Whole-graph log score (sum over nodes), with a per-node cache
graph_logscore <- function(DAG, tXX, n, a, U, approx = FALSE, cache = NULL) {
  q <- ncol(DAG)
  if (is.null(cache)) cache <- rep(NA_real_, q)
  for (j in 1:q) if (is.na(cache[j]))
    cache[j] <- node_score(j, DAG, tXX, n, a, U, approx)
  list(total = sum(cache), cache = cache)
}

## ----------------------------------------------------------------------------
## 2.  Neighbourhood operators (Insert / Delete / Reverse)
## ----------------------------------------------------------------------------
apply_op <- function(op, A, x, y) {
  if (op == 1) A[x, y] <- 1
  if (op == 2) A[x, y] <- 0
  if (op == 3) { A[x, y] <- 0; A[y, x] <- 1 }
  A
}

## Enumerate all valid single-edge moves; returns a matrix with columns
## (op, x, y) and, for each, the set of affected child nodes.
enumerate_moves <- function(DAG) {
  q <- ncol(DAG); A_na <- DAG; diag(A_na) <- NA
  ins <- which(A_na == 0, TRUE); ins <- if (length(ins)) cbind(1, ins) else NULL
  del <- which(A_na == 1, TRUE); del <- if (length(del)) cbind(2, del) else NULL
  rev <- which(A_na == 1, TRUE); rev <- if (length(rev)) cbind(3, rev) else NULL
  O <- rbind(ins, del, rev)
  keep <- logical(nrow(O))
  for (i in seq_len(nrow(O)))
    keep[i] <- gRbase::is.DAG(apply_op(O[i,1], DAG, O[i,2], O[i,3]))
  O[keep, , drop = FALSE]
}

## Affected child node(s) of a move
affected_nodes <- function(op, x, y) if (op == 3) c(x, y) else y

## ----------------------------------------------------------------------------
## 3.  Random-walk Metropolis-Hastings (baseline)
## ----------------------------------------------------------------------------
sampler_randomwalk <- function(S, burn, data, a, U, w, approx = FALSE) {
  X <- scale(data, scale = FALSE); tXX <- crossprod(X)
  n <- nrow(data); q <- ncol(data)
  DAG <- matrix(0, q, q)
  gs  <- graph_logscore(DAG, tXX, n, a, U, approx)
  cur_cache <- gs$cache; cur_score <- gs$total
  size <- numeric(S + burn)
  Graphs <- array(0, c(q, q, S))
  logprior <- function(op) c(log(w/(1-w)), log((1-w)/w), 0)[op]
  for (it in 1:(S + burn)) {
    O <- enumerate_moves(DAG); m <- nrow(O)
    i <- sample(m, 1); op <- O[i,1]; x <- O[i,2]; y <- O[i,3]
    prop <- apply_op(op, DAG, x, y)
    aff  <- affected_nodes(op, x, y)
    prop_cache <- cur_cache
    for (nd in aff) prop_cache[nd] <- node_score(nd, prop, tXX, n, a, U, approx)
    prop_score <- sum(prop_cache)
    Oprop <- enumerate_moves(prop)
    log_acc <- (prop_score - cur_score) + logprior(op) +
               (log(m) - log(nrow(Oprop)))
    if (log(runif(1)) < log_acc) {
      DAG <- prop; cur_cache <- prop_cache; cur_score <- prop_score
    }
    size[it] <- sum(DAG)
    if (it > burn) Graphs[,,it-burn] <- DAG
  }
  list(Graphs = Graphs, size = size[(burn+1):(burn+S)], size_full = size)
}

## ----------------------------------------------------------------------------
## 4.  Locally-balanced informed sampler (Algorithm 1)
## ----------------------------------------------------------------------------
##  Balancing functions g(t) with g(t) = t g(1/t):
##    sqrt    : g(t) = sqrt(t)              (locally balanced)
##    barker  : g(t) = t / (1 + t)
balancing <- function(logr, type = "sqrt") {
  if (type == "sqrt")   return(exp(0.5 * logr))
  if (type == "barker") return(1 / (1 + exp(-logr)))
  stop("unknown balancing function")
}

## Compute, for every neighbour of DAG, the log score ratio relative to DAG,
## reusing the current node-score cache.  Returns the move table O, the vector
## of log ratios, and Z_g (the neighbourhood normaliser).
neighbour_weights <- function(DAG, tXX, n, a, U, w, cur_cache, cur_score,
                              g = "sqrt", approx = FALSE) {
  O <- enumerate_moves(DAG); m <- nrow(O)
  logr <- numeric(m)
  logprior <- function(op) c(log(w/(1-w)), log((1-w)/w), 0)[op]
  for (i in 1:m) {
    op <- O[i,1]; x <- O[i,2]; y <- O[i,3]
    prop <- apply_op(op, DAG, x, y)
    aff  <- affected_nodes(op, x, y)
    delta <- 0
    for (nd in aff)
      delta <- delta + (node_score(nd, prop, tXX, n, a, U, approx) - cur_cache[nd])
    logr[i] <- delta + logprior(op)        # log pi(D')/pi(D)
  }
  gvals <- balancing(logr, g)
  list(O = O, logr = logr, gvals = gvals, Zg = sum(gvals))
}

sampler_informed <- function(S, burn, data, a, U, w,
                             g = "sqrt", approx = FALSE) {
  X <- scale(data, scale = FALSE); tXX <- crossprod(X)
  n <- nrow(data); q <- ncol(data)
  DAG <- matrix(0, q, q)
  gs  <- graph_logscore(DAG, tXX, n, a, U, approx)
  cur_cache <- gs$cache; cur_score <- gs$total
  nb <- neighbour_weights(DAG, tXX, n, a, U, w, cur_cache, cur_score, g, approx)
  size <- numeric(S + burn); Graphs <- array(0, c(q, q, S))

  for (it in 1:(S + burn)) {
    ## propose a neighbour with probability g_i / Z_g
    probs <- nb$gvals / nb$Zg
    i <- sample(length(probs), 1, prob = probs)
    op <- nb$O[i,1]; x <- nb$O[i,2]; y <- nb$O[i,3]
    prop <- apply_op(op, DAG, x, y)
    aff  <- affected_nodes(op, x, y)
    prop_cache <- cur_cache
    for (nd in aff) prop_cache[nd] <- node_score(nd, prop, tXX, n, a, U, approx)
    prop_score <- sum(prop_cache)

    ## neighbourhood normaliser at the proposed graph
    nb_prop <- neighbour_weights(prop, tXX, n, a, U, w, prop_cache, prop_score,
                                 g, approx)
    ## acceptance = min(1, Z_g(D) / Z_g(D'))   (Proposition 3)
    log_acc <- log(nb$Zg) - log(nb_prop$Zg)
    if (log(runif(1)) < log_acc) {
      DAG <- prop; cur_cache <- prop_cache; cur_score <- prop_score; nb <- nb_prop
    }
    size[it] <- sum(DAG)
    if (it > burn) Graphs[,,it-burn] <- DAG
  }
  list(Graphs = Graphs, size = size[(burn+1):(burn+S)], size_full = size)
}

## ----------------------------------------------------------------------------
## 5.  Posterior summaries and efficiency metrics
## ----------------------------------------------------------------------------
get_edgeprobs <- function(out) apply(out$Graphs, c(1,2), mean)
get_MPMdag    <- function(out) (get_edgeprobs(out) > 0.5) * 1

ess_per_sec <- function(out, seconds) {
  e <- as.numeric(coda::effectiveSize(coda::mcmc(out$size)))
  c(ESS = e, ESS_per_sec = e / seconds)
}

## ----------------------------------------------------------------------------
## 6.  Example driver: compare random-walk vs informed on one simulated DAG
## ----------------------------------------------------------------------------
##  Uncomment to run.
##
##  set.seed(1)
##  q <- 30; n <- 300; w <- 0.15
##  L1 <- huge::huge.generator(d = q, prob = w, vis = FALSE)
##  DAGtrue <- as.matrix(igraph::as_adjacency_matrix(
##               igraph::graph_from_adjacency_matrix(L1$theta)))
##  L <- matrix(runif(q*q, 1, 2) * sample(c(-1,1), q*q, TRUE), q, q) * DAGtrue
##  diag(L) <- 1; Omega <- t(L) %*% L
##  X <- mvtnorm::rmvnorm(n, sigma = solve(Omega))
##
##  t_rw  <- system.time(
##    out_rw  <- sampler_randomwalk(S = 20000, burn = 5000, data = X,
##                                  a = q, U = diag(q), w = w))[3]
##  t_inf <- system.time(
##    out_inf <- sampler_informed(S = 20000, burn = 5000, data = X,
##                                a = q, U = diag(q), w = w,
##                                g = "sqrt", approx = TRUE))[3]
##
##  cat("Random-walk:", ess_per_sec(out_rw,  t_rw),  "\n")
##  cat("Informed   :", ess_per_sec(out_inf, t_inf), "\n")

cat("nCPNG informed-sampler module loaded. See driver block at end of file.\n")
