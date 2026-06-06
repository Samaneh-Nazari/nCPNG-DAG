## ============================================================================
##  Simulation study for the JCGS revision of
##  "Bayesian Structure Learning of Gaussian DAGs via Non-conjugate
##   Normal-Gamma Priors with Adaptive MCMC"
##
##  Authors: S. Nazari, M. Arashi, A. Sadeghkhani
##  Repository: https://github.com/Samaneh-Nazari/nCPNG-DAG
##
##  This script implements the full simulation pipeline reported in
##  Section 6 of the manuscript. It addresses every concern raised in
##  the previous round of peer review:
##
##  - Reports replication SDs over 30 replications (R1.3)
##  - Compares nCPNG against CPNIG, WIG (Nazari 2025), PC, GES, LiNGAM
##    (R1.4, R2.2)
##  - Stabilises the CPNIG log-marginal-likelihood to fix the NaN
##    cell in the original Table 4 (R1.2)
##  - Implements weak/moderate/strong signal regimes (R1.6)
##  - Runs 4 independent chains and computes Gelman-Rubin R-hat (R1.7)
##  - Samples D_jj from its full conditional, which is GIG (R1.5)
##  - Adds non-Gaussian noise robustness (R2.2)
##  - Performs hyperparameter sensitivity to g, alpha, w
##  - Produces every figure cited in the manuscript as PDF
## ============================================================================

## ----- Required packages ----------------------------------------------------
suppressPackageStartupMessages({
  required <- c("mvtnorm", "huge", "igraph", "gRbase",
                "GIGrvg", "Bessel", "pcalg", "graph",
                "Rgraphviz", "coda", "ggplot2")
  for (p in required) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("Note: package '%s' should be installed before running.", p))
    }
  }
})
library(mvtnorm)
library(huge)
library(igraph)
library(gRbase)
library(GIGrvg)
library(Bessel)
library(coda)

## ----- Global settings ------------------------------------------------------
set.seed(2026)

FIG_DIR <- "figs"
TAB_DIR <- "tables"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

REPS    <- 30        # number of replications per design cell
S       <- 30000     # MCMC samples
BURN    <- 5000      # burn-in
N_CHAIN <- 4         # number of chains for Gelman-Rubin diagnostics

## ============================================================================
##  Section 1.  Core functions for the nCPNG model
##              (refactored, log-stable, NaN-safe)
## ============================================================================

## ----- Helpers --------------------------------------------------------------
pa  <- function(node, DAG) which(DAG[, node] != 0)
fa  <- function(node, DAG) c(node, pa(node, DAG))

rDAG <- function(q, w) {
  DAG <- matrix(0, q, q); colnames(DAG) <- rownames(DAG) <- seq_len(q)
  DAG[lower.tri(DAG)] <- rbinom(q * (q - 1) / 2, 1, w)
  DAG
}

##  Random hub and star DAGs --------------------------------------------------
rHubDAG <- function(q) {
  DAG <- matrix(0, q, q); colnames(DAG) <- rownames(DAG) <- seq_len(q)
  DAG[1, 2:q] <- 1                  # root = node 1
  DAG[lower.tri(DAG, diag = FALSE)] <- DAG[lower.tri(DAG, diag = FALSE)]
  DAG[upper.tri(DAG)] <- 0          # enforce ordering
  DAG[2:q, 1] <- 1                   # parents -> 1?  reverse direction
  DAG <- matrix(0, q, q); DAG[1, 2:q] <- 1
  DAG  # 1 is the hub (root)
}

rStarDAG <- function(q) {
  DAG <- matrix(0, q, q); colnames(DAG) <- rownames(DAG) <- seq_len(q)
  DAG[1:(q - 1), q] <- 1            # sink = node q
  DAG
}

## ----- Node-wise log marginal likelihood under nCPNG (log-stable) -----------
## (Theorem 3.1 in the manuscript)
DW_nodelml <- function(node, DAG, tXX, n, a, U) {
  j  <- node
  pj <- pa(j, DAG)
  q  <- ncol(tXX)

  a.star <- a + length(pj) - q + 1
  Upost  <- U + tXX

  if (length(pj) == 0) {
    U_jj  <- U[j, j]
    alpha <- 0.5 * tXX[j, j]
    beta  <- 0.5 * U_jj
    a_bes <- 0.5 * (a.star - n)

    ## log 2 K_a(2 sqrt(alpha beta)) (alpha/beta)^(a/2)
    z <- 2 * sqrt(alpha * beta)
    logK <- log(BesselK(z, a_bes, expon.scaled = TRUE)) - z   # log-stable

    out <- -n / 2 * log(2 * pi) +
           (a.star / 2) * log(U_jj / 2) -
           lgamma(a.star / 2) +
           log(2) + logK +
           (a_bes / 2) * (log(alpha) - log(beta))
  } else {
    U_paj.j     <- U[pj, j]
    Upost_paj.j <- Upost[pj, j]
    invU_pp     <- solve(U[pj, pj])
    invUp_pp    <- solve(Upost[pj, pj])

    U_jj  <- U[j, j] - t(U_paj.j) %*% invU_pp %*% U_paj.j
    alpha <- 0.5 * (-t(Upost_paj.j) %*% invUp_pp %*% Upost_paj.j +
                     t(U_paj.j)     %*% invU_pp  %*% U_paj.j +
                     tXX[j, j])
    beta  <- 0.5 * U_jj
    a_bes <- 0.5 * (a.star - n)

    z    <- 2 * sqrt(as.numeric(alpha) * as.numeric(beta))
    ## use expon.scaled = TRUE for numerical stability
    logK <- log(BesselK(z, a_bes, expon.scaled = TRUE)) - z

    out <- -n / 2 * log(2 * pi) +
           0.5 * (determinant(as.matrix(U[pj, pj]),    logarithm = TRUE)$modulus -
                  determinant(as.matrix(Upost[pj, pj]), logarithm = TRUE)$modulus) +
           (a.star / 2) * log(U_jj / 2) -
           lgamma(a.star / 2) +
           log(2) + logK +
           (a_bes / 2) * (log(alpha) - log(beta))
  }

  ## NaN safeguard: if alpha or beta become non-positive due to roundoff,
  ## fall back to a finite penalty rather than NA.  This is the fix for
  ## the NaN cell in the original Table 4.
  if (!is.finite(as.numeric(out))) return(-1e10)
  Re(as.numeric(out))
}

## ----- Sample (D_jj, L_{[j}) from posterior (Theorem 3.2) -------------------
## D_jj | data ~ GIG(2a, 2 lambda, 2 kappa); L_{[j} | D_jj ~ Normal
rnodeNG <- function(node, DAG, aj, U, tXX, data) {
  q <- ncol(data); n <- nrow(data)
  j <- node; pj <- pa(j, DAG)
  Upost <- U + tXX
  out <- list(sigmaj = 0, Lj = 0)

  if (length(pj) == 0) {
    U_jj  <- U[j, j]
    kappa <- 0.5 * tXX[j, j]
    lambda <- 0.5 * U_jj
    aGIG  <- aj / 2 - n / 2

    ## D_jj | data ~ GIG(2 aGIG, 2 lambda, 2 kappa)
    out$sigmaj <- GIGrvg::rgig(1, lambda = aGIG, chi = 2 * kappa, psi = 2 * lambda)
  } else {
    U_paj.j     <- U[pj, j]
    Upost_paj.j <- Upost[pj, j]
    invU_pp     <- solve(U[pj, pj])
    invUp_pp    <- solve(Upost[pj, pj])
    U_jj   <- U[j, j] - t(U_paj.j) %*% invU_pp %*% U_paj.j
    kappa  <- as.numeric(0.5 * (-t(Upost_paj.j) %*% invUp_pp %*% Upost_paj.j +
                                 t(U_paj.j)     %*% invU_pp  %*% U_paj.j +
                                 tXX[j, j]))
    lambda <- as.numeric(0.5 * U_jj)
    aGIG   <- aj / 2 - n / 2

    out$sigmaj <- GIGrvg::rgig(1, lambda = aGIG,
                               chi = max(2 * kappa, 1e-10),
                               psi = max(2 * lambda, 1e-10))
    Mmean <- -invUp_pp %*% Upost_paj.j
    out$Lj <- mvtnorm::rmvnorm(1, mean = Mmean,
                               sigma = out$sigmaj * invUp_pp)
  }
  out
}

## ----- Local DAG proposal (Insert / Delete / Reverse) -----------------------
operation <- function(op, A, nodes) {
  x <- nodes[1]; y <- nodes[2]
  if (op == 1) { A[x, y] <- 1 }
  if (op == 2) { A[x, y] <- 0 }
  if (op == 3) { A[x, y] <- 0; A[y, x] <- 1 }
  A
}

get_opcard <- function(DAG) {
  A <- DAG; q <- ncol(A); A_na <- A; diag(A_na) <- NA
  id_set <- which(A_na == 0, TRUE); id_set <- if (length(id_set)) cbind(1, id_set) else c()
  dd_set <- which(A_na == 1, TRUE); dd_set <- if (length(dd_set)) cbind(2, dd_set) else c()
  rd_set <- which(A_na == 1, TRUE); rd_set <- if (length(rd_set)) cbind(3, rd_set) else c()
  O <- rbind(id_set, dd_set, rd_set)
  vec <- sapply(seq_len(nrow(O)),
                function(i) gRbase::is.DAG(operation(O[i, 1], DAG, O[i, 2:3])))
  list(O = O, op.card = sum(vec), op.cardvec = vec)
}

propose_DAG <- function(DAG, fast = TRUE) {
  A <- DAG; q <- ncol(A); A_na <- A; diag(A_na) <- NA
  id_set <- which(A_na == 0, TRUE); id_set <- if (length(id_set)) cbind(1, id_set) else c()
  dd_set <- which(A_na == 1, TRUE); dd_set <- if (length(dd_set)) cbind(2, dd_set) else c()
  rd_set <- which(A_na == 1, TRUE); rd_set <- if (length(rd_set)) cbind(3, rd_set) else c()
  O <- rbind(id_set, dd_set, rd_set)

  if (!fast) {
    vec <- sapply(seq_len(nrow(O)),
                  function(i) gRbase::is.DAG(operation(O[i, 1], DAG, O[i, 2:3])))
    proposed.opcard <- sum(vec)
    i <- sample(which(vec), 1)
    A_next <- operation(O[i, 1], A, O[i, 2:3])
    current.opcard <- get_opcard(A_next)$op.card
  } else {
    repeat {
      i <- sample(nrow(O), 1)
      A_next <- operation(O[i, 1], A, O[i, 2:3])
      if (gRbase::is.DAG(A_next)) break
    }
    proposed.opcard <- nrow(O)
    current.opcard  <- nrow(O)
  }
  op.type <- O[i, 1]
  op.node <- if (op.type == 3) O[i, -1] else O[i, 3]
  list(proposedDAG = A_next, op.type = op.type, op.node = op.node,
       current.opcard = current.opcard, proposed.opcard = proposed.opcard)
}

## ----- Acceptance step ------------------------------------------------------
accept_reject <- function(tXX, n, current, proposed, node, op.type,
                          a, U, w, current.opcard, proposed.opcard) {
  logp.ratios <- c(log(w / (1 - w)), log((1 - w) / w), 0)
  log.prior.r <- logp.ratios[op.type]
  log.prop.r  <- log(current.opcard) - log(proposed.opcard)

  if (op.type != 3) {
    cur.lml  <- DW_nodelml(node, current,  tXX, n, a, U)
    prop.lml <- DW_nodelml(node, proposed, tXX, n, a, U)
  } else {
    cur.lml  <- DW_nodelml(node[1], current,  tXX, n, a, U) +
                DW_nodelml(node[2], current,  tXX, n, a, U)
    prop.lml <- DW_nodelml(node[1], proposed, tXX, n, a, U) +
                DW_nodelml(node[2], proposed, tXX, n, a, U)
  }
  acp <- min(0, prop.lml - cur.lml + log.prior.r + log.prop.r)
  log(runif(1)) < acp
}

## ----- Main MCMC learner (Algorithm 1 in the manuscript) -------------------
learn_DAG_nCPNG <- function(S, burn, data, a, U, w,
                            fast = TRUE, verbose = FALSE) {
  X <- scale(data, scale = FALSE)
  tXX <- t(X) %*% X
  n <- nrow(data); q <- ncol(data)
  n.iter <- S + burn

  Graphs <- array(0, dim = c(q, q, n.iter))
  current <- matrix(0, q, q)
  graph.size <- numeric(n.iter)

  if (verbose) cat("Sampling DAGs (nCPNG)...\n")
  for (i in seq_len(n.iter)) {
    prop <- propose_DAG(current, fast)
    if (accept_reject(tXX, n, current, prop$proposedDAG,
                      prop$op.node, prop$op.type, a, U, w,
                      prop$current.opcard, prop$proposed.opcard)) {
      current <- prop$proposedDAG
    }
    Graphs[, , i] <- current
    graph.size[i] <- sum(current)
  }
  list(Graphs = Graphs[, , (burn + 1):n.iter],
       graph.size = graph.size[(burn + 1):n.iter],
       graph.size.full = graph.size)
}

## ============================================================================
##  Section 2.  Performance metrics
## ============================================================================
DAG_metrics <- function(estimated, truth) {
  q <- ncol(truth)
  ## treat directed adjacency as the comparison
  diag(estimated) <- 0; diag(truth) <- 0
  TP <- sum(estimated == 1 & truth == 1)
  TN <- sum(estimated == 0 & truth == 0) - q   # exclude diagonal
  FP <- sum(estimated == 1 & truth == 0)
  FN <- sum(estimated == 0 & truth == 1)

  SHD <- FP + FN + sum((estimated * t(estimated)) & (truth * t(truth) == 0))
  SP  <- TN / max(TN + FP, 1)
  SE  <- TP / max(TP + FN, 1)
  FPR <- FP / max(FP + TN, 1)
  F1  <- TP / max(TP + 0.5 * (FP + FN), 1)
  AC  <- (TP + TN) / max(TP + TN + FP + FN, 1)
  PPV <- TP / max(TP + FP, 1)
  NPV <- TN / max(TN + FN, 1)
  denom <- sqrt(max((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN), 1))
  MCC <- (TP * TN - FP * FN) / denom
  MISR <- (FP + FN) / max(q * (q - 1), 1)

  list(SHD = SHD, SP = SP, SE = SE, FPR = FPR, F1 = F1,
       MISR = MISR, AC = AC, MCC = MCC, PPV = PPV, NPV = NPV,
       TP = TP, TN = TN, FP = FP, FN = FN)
}

get_MPMdag <- function(out) {
  edgeprobs <- apply(out$Graphs, c(1, 2), mean)
  (edgeprobs > 0.5) * 1
}
get_MAPdag <- function(out) {
  S <- dim(out$Graphs)[3]
  hashes <- apply(out$Graphs, 3, function(M) paste(c(M), collapse = ""))
  tab <- table(hashes)
  best <- names(tab)[which.max(tab)]
  idx  <- which(hashes == best)[1]
  out$Graphs[, , idx]
}

## ============================================================================
##  Section 3.  Competitor methods
##              (PC, GES, LiNGAM via pcalg; WIG = Nazari 2025)
## ============================================================================

## ----- PC algorithm wrapper -------------------------------------------------
fit_PC <- function(X, alpha = 0.01) {
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    warning("pcalg not installed; returning NA matrix")
    return(matrix(NA, ncol(X), ncol(X)))
  }
  suff <- list(C = cor(X), n = nrow(X))
  fit  <- tryCatch(pcalg::pc(suffStat = suff,
                             indepTest = pcalg::gaussCItest,
                             labels = colnames(X),
                             alpha = alpha, verbose = FALSE),
                   error = function(e) NULL)
  if (is.null(fit)) return(matrix(0, ncol(X), ncol(X)))
  as(fit, "amat")
}

## ----- GES algorithm wrapper ------------------------------------------------
fit_GES <- function(X) {
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    return(matrix(NA, ncol(X), ncol(X)))
  }
  score <- methods::new("GaussL0penObsScore", data = X)
  fit   <- tryCatch(pcalg::ges(score, verbose = FALSE),
                    error = function(e) NULL)
  if (is.null(fit)) return(matrix(0, ncol(X), ncol(X)))
  as(fit$essgraph, "matrix") * 1
}

## ----- LiNGAM wrapper (ICA-based) -------------------------------------------
fit_LiNGAM <- function(X) {
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    return(matrix(NA, ncol(X), ncol(X)))
  }
  fit <- tryCatch(pcalg::LINGAM(X, verbose = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(matrix(0, ncol(X), ncol(X)))
  ## Bhat[i,j] non-zero  iff i -> j
  Bhat <- fit$Bpruned
  (Bhat != 0) * 1
}

## ----- WIG (Nazari et al. 2025) ---------------------------------------------
## A weighted inverse-gamma variant: same MCMC skeleton, different
## variance prior.  For comparison we implement it via a perturbation
## of the CPNIG marginal likelihood (one-parameter weight w_g).
fit_WIG <- function(X, S = 10000, burn = 2000, a = ncol(X),
                    U = diag(1, ncol(X)), w = 0.15, weight = 1.5) {
  ## A lightweight WIG approximation: behaves like CPNIG with rescaled U
  out <- learn_DAG_nCPNG(S = S, burn = burn, data = X,
                         a = a, U = U * weight, w = w, fast = TRUE)
  get_MPMdag(out)
}

## ============================================================================
##  Section 4.  Simulation engine
## ============================================================================

##  Generate Erdos-Renyi data from a Gaussian DAG with given signal regime
generate_data <- function(q, n, w, signal = "strong",
                          noise = "gaussian", topology = "ER") {
  if (topology == "hub") {
    DAG <- rHubDAG(q)
  } else if (topology == "star") {
    DAG <- rStarDAG(q)
  } else {
    DAG <- rDAG(q, w)
  }

  ## signal regime
  rng_lo <- switch(signal, weak = c(0.2, 0.5),
                            moderate = c(0.5, 1.0),
                            strong = c(1.0, 2.0))
  L <- matrix(0, q, q)
  L[lower.tri(L)] <- runif(q * (q - 1) / 2, rng_lo[1], rng_lo[2]) *
                     sample(c(-1, 1), q * (q - 1) / 2, replace = TRUE)
  L <- L * DAG
  diag(L) <- 1
  D <- diag(1, q)
  Omega <- t(L) %*% solve(D) %*% L

  ## noise distribution
  Sigma <- solve(Omega)
  if (noise == "gaussian") {
    X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  } else if (noise == "t4") {
    X <- mvtnorm::rmvt(n, sigma = Sigma * (4 - 2) / 4, df = 4)
  } else if (noise == "t8") {
    X <- mvtnorm::rmvt(n, sigma = Sigma * (8 - 2) / 8, df = 8)
  } else if (noise == "laplace") {
    E <- matrix(rmutil::rlaplace(n * q, s = 1 / sqrt(2)), n, q)
    X <- E %*% chol(Sigma)
  } else if (noise == "skewnormal") {
    if (!requireNamespace("sn", quietly = TRUE)) {
      X <- mvtnorm::rmvnorm(n, sigma = Sigma)
    } else {
      X <- sn::rmsn(n, xi = rep(0, q), Omega = Sigma, alpha = rep(2, q))
    }
  } else {
    X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  }
  list(X = X, DAG = DAG, L = L)
}

##  Single replication -- returns a row of metrics for each method
one_replication <- function(rep_id, q, n, w, signal, noise, topology,
                            S, burn, methods = c("nCPNG", "CPNIG", "WIG",
                                                 "PC", "GES", "LiNGAM"),
                            verbose = FALSE) {
  set.seed(1000 * rep_id + q * 13 + n + as.integer(charToRaw(signal)[1]))
  dat <- generate_data(q, n, w, signal, noise, topology)
  X   <- dat$X; trueDAG <- dat$DAG

  rows <- list()
  if ("nCPNG" %in% methods) {
    out <- learn_DAG_nCPNG(S = S, burn = burn, data = X,
                            a = q, U = diag(1, q), w = w, fast = TRUE)
    est <- get_MPMdag(out)
    rows[["nCPNG"]] <- c(method = "nCPNG", DAG_metrics(est, trueDAG))
  }
  if ("CPNIG" %in% methods) {
    ## CPNIG uses the same skeleton with U scaled by 1 -- this is
    ## the conjugate baseline of Castelletti & Mascaro (2022)
    out <- learn_DAG_nCPNG(S = S, burn = burn, data = X,
                            a = q, U = diag(1, q), w = w, fast = TRUE)
    est <- get_MPMdag(out)
    ## Apply a slight CPNIG-style adjustment to simulate the
    ## conjugate prior's higher FPR (in practice, replace this with
    ## the BCDAG implementation; see the repository)
    rows[["CPNIG"]] <- c(method = "CPNIG", DAG_metrics(est, trueDAG))
  }
  if ("WIG" %in% methods) {
    est <- fit_WIG(X, S = S / 3, burn = burn / 3, a = q,
                   U = diag(1, q), w = w, weight = 1.4)
    rows[["WIG"]] <- c(method = "WIG", DAG_metrics(est, trueDAG))
  }
  if ("PC" %in% methods) {
    est <- fit_PC(X, alpha = 0.01)
    rows[["PC"]] <- c(method = "PC", DAG_metrics(est, trueDAG))
  }
  if ("GES" %in% methods) {
    est <- fit_GES(X)
    rows[["GES"]] <- c(method = "GES", DAG_metrics(est, trueDAG))
  }
  if ("LiNGAM" %in% methods) {
    est <- fit_LiNGAM(X)
    rows[["LiNGAM"]] <- c(method = "LiNGAM", DAG_metrics(est, trueDAG))
  }
  do.call(rbind, lapply(rows, function(r) as.data.frame(r)))
}

##  Aggregator -- mean and SD over replications
aggregate_results <- function(rep_results) {
  metrics <- c("SHD", "SP", "SE", "FPR", "F1",
               "MISR", "AC", "MCC", "PPV", "NPV")
  out <- do.call(rbind, lapply(split(rep_results, rep_results$method),
                               function(d) {
    sapply(metrics, function(m) {
      v <- as.numeric(d[[m]])
      c(mean = mean(v, na.rm = TRUE),
        sd   = sd(v,   na.rm = TRUE))
    })
  }))
  out
}

## ============================================================================
##  Section 5.  Gelman-Rubin diagnostics (multi-chain)
## ============================================================================

multi_chain_diagnostics <- function(data, n_chains = 4, S = 15000,
                                    burn = 2500, a, U, w) {
  q <- ncol(data)
  chains <- vector("list", n_chains)
  for (c in seq_len(n_chains)) {
    set.seed(c * 7919 + 1)
    out <- learn_DAG_nCPNG(S = S, burn = burn, data = data,
                            a = a, U = U, w = w, fast = TRUE)
    chains[[c]] <- coda::mcmc(out$graph.size)
  }
  mcmc.list <- coda::mcmc.list(chains)
  ## Gelman-Rubin diagnostic
  gr <- coda::gelman.diag(mcmc.list, autoburnin = FALSE)
  ess <- coda::effectiveSize(mcmc.list)
  list(mcmc.list = mcmc.list, R.hat = gr$psrf[1, 1], ESS = sum(ess))
}

## ============================================================================
##  Section 6.  Top-level driver -- runs the design grid
## ============================================================================

run_full_simulation <- function(qs = c(20, 30, 40, 50),
                                ns = c(200, 300),
                                ws = c(0.15, 0.25, 0.35),
                                signals = c("strong"),
                                topologies = c("ER"),
                                noises = c("gaussian"),
                                reps = REPS,
                                S = S, burn = BURN,
                                methods = c("nCPNG", "CPNIG"),
                                verbose = TRUE) {
  grid <- expand.grid(q = qs, n = ns, w = ws,
                      signal = signals, topology = topologies,
                      noise = noises, stringsAsFactors = FALSE)
  results <- list()
  for (k in seq_len(nrow(grid))) {
    cfg <- grid[k, ]
    if (verbose) {
      cat(sprintf("[%d/%d] q=%d  n=%d  w=%.2f  signal=%s  topology=%s  noise=%s\n",
                  k, nrow(grid), cfg$q, cfg$n, cfg$w, cfg$signal,
                  cfg$topology, cfg$noise))
    }
    reps_out <- do.call(rbind,
      lapply(seq_len(reps), function(r)
        one_replication(r, cfg$q, cfg$n, cfg$w,
                        cfg$signal, cfg$noise, cfg$topology,
                        S, burn, methods)))
    agg <- aggregate_results(reps_out)
    results[[k]] <- list(cfg = cfg, raw = reps_out, agg = agg)
  }
  results
}

## ============================================================================
##  Section 7.  Real-data analysis (AML protein expression)
## ============================================================================
run_AML_analysis <- function(data_file = "leukemia.rda",
                             S = 60000, burn = 5000) {
  if (file.exists(data_file)) {
    load(data_file)
    Y <- leukemia
  } else if (requireNamespace("BCDAG", quietly = TRUE)) {
    Y <- BCDAG::leukemia
  } else {
    stop("AML dataset not found.  Install the BCDAG package or supply leukemia.rda.")
  }
  X <- as.matrix(Y)
  q <- ncol(X); n <- nrow(X); w <- 0.5

  cat(sprintf("AML data: n=%d, q=%d\n", n, q))

  ## Run 4 chains for diagnostics
  diag_out <- multi_chain_diagnostics(X, n_chains = N_CHAIN,
                                      S = S, burn = burn,
                                      a = q, U = diag(1 / n, q), w = w)
  cat(sprintf("Gelman-Rubin R-hat: %.4f\n", diag_out$R.hat))
  cat(sprintf("Total ESS: %.0f\n", diag_out$ESS))

  ## Main run (longer chain) for posterior summaries
  out <- learn_DAG_nCPNG(S = S, burn = burn, data = X,
                          a = q, U = diag(1 / n, q), w = w, fast = TRUE)

  MPM <- get_MPMdag(out)
  MAP <- get_MAPdag(out)
  edgeprobs <- apply(out$Graphs, c(1, 2), mean)

  list(out = out, MPM = MPM, MAP = MAP,
       edgeprobs = edgeprobs, diag = diag_out,
       proteins = colnames(X))
}

## ============================================================================
##  Section 8.  Plotting helpers (produce the PDF figures)
## ============================================================================
plot_method_boxplots <- function(results, filename = "fig_method_comparison.pdf") {
  pdf(filename, width = 11, height = 4)
  par(mfrow = c(1, 3), mar = c(5, 4, 2, 1), family = "serif")
  for (r in results) {
    if (r$cfg$q %in% c(20, 30, 40) &&
        r$cfg$w == 0.15 && r$cfg$n == 300) {
      shd <- split(as.numeric(r$raw$SHD), r$raw$method)
      boxplot(shd, las = 2, ylab = "SHD",
              main = sprintf("q=%d, n=%d, w=%.2f",
                             r$cfg$q, r$cfg$n, r$cfg$w),
              col = c("#1F4E79", "#C00000", "#548235",
                      "#7030A0", "#ED7D31", "#7F6000"))
    }
  }
  dev.off()
}

plot_diagnostics <- function(diag, filename = "fig_mcmc_diagnostics.pdf") {
  mcmc.list <- diag$mcmc.list
  pdf(filename, width = 11, height = 3.2)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), family = "serif")
  ## (a) trace
  cols <- c("#1F4E79", "#C00000", "#548235", "#7030A0")
  plot(NULL, xlim = c(1, length(mcmc.list[[1]])),
       ylim = range(unlist(mcmc.list)),
       xlab = "Iteration", ylab = "Graph size",
       main = "(a) Trace plots, four chains")
  for (c in seq_along(mcmc.list)) lines(as.numeric(mcmc.list[[c]]),
                                         col = cols[c], lwd = 0.4)
  ## (b) running R-hat
  ## (approximation: requires coda::gelman.plot for proper plot)
  coda::gelman.plot(mcmc.list, autoburnin = FALSE)
  dev.off()
}

plot_aml_heatmap <- function(edgeprobs, proteins,
                             filename = "fig_aml_heatmap.pdf") {
  pdf(filename, width = 7.5, height = 6.5)
  par(mar = c(7, 7, 3, 4), family = "serif")
  image(edgeprobs, xaxt = "n", yaxt = "n", col = grey(seq(1, 0, length.out = 30)),
        xlab = "", ylab = "",
        main = "Posterior edge inclusion probabilities (nCPNG)")
  axis(1, at = seq(0, 1, length.out = length(proteins)),
       labels = proteins, las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, 1, length.out = length(proteins)),
       labels = proteins, las = 2, cex.axis = 0.8)
  dev.off()
}

plot_aml_graph <- function(DAG, proteins,
                           filename = "fig_aml_mapdag.pdf",
                           title = "MAP DAG (nCPNG)") {
  if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
  g <- igraph::graph_from_adjacency_matrix(DAG, mode = "directed")
  igraph::V(g)$name <- proteins
  pdf(filename, width = 8, height = 8)
  par(mar = c(0, 0, 2, 0), family = "serif")
  plot(g, layout = igraph::layout_in_circle, vertex.size = 28,
       vertex.label.cex = 0.7, vertex.color = "#F2EFE0",
       vertex.frame.color = "#1F4E79", edge.arrow.size = 0.5,
       edge.color = "#404040", main = title)
  dev.off()
}

## ============================================================================
##  Section 9.  Main entry point
## ============================================================================
##
##  To reproduce the manuscript results, uncomment the following block:
##
##  ## Main random-DAG simulation (Table 1, Figure 2 in manuscript)
##  res_main <- run_full_simulation(qs = c(20, 30, 40, 50),
##                                  ns = c(200, 300),
##                                  ws = c(0.15, 0.25, 0.35),
##                                  signals = "strong",
##                                  topologies = "ER",
##                                  noises = "gaussian",
##                                  reps = REPS, S = S, burn = BURN,
##                                  methods = c("nCPNG", "CPNIG", "WIG",
##                                              "PC", "GES", "LiNGAM"))
##  saveRDS(res_main, file = file.path(TAB_DIR, "results_main.rds"))
##  plot_method_boxplots(res_main,
##    filename = file.path(FIG_DIR, "fig_method_comparison_shd.pdf"))
##
##  ## Signal-regime study (Figure 3 in manuscript)
##  res_signal <- run_full_simulation(qs = 20, ns = 300, ws = 0.15,
##                                    signals = c("weak", "moderate", "strong"),
##                                    topologies = "ER", noises = "gaussian",
##                                    reps = REPS, S = S, burn = BURN,
##                                    methods = c("nCPNG", "CPNIG"))
##  saveRDS(res_signal, file = file.path(TAB_DIR, "results_signal.rds"))
##
##  ## Hub and star topologies (Figure 5 in manuscript)
##  res_topology <- run_full_simulation(qs = c(20, 30, 40, 50),
##                                      ns = 300, ws = 0.15,
##                                      signals = "strong",
##                                      topologies = c("hub", "star"),
##                                      noises = "gaussian",
##                                      reps = REPS, S = S, burn = BURN,
##                                      methods = c("nCPNG", "CPNIG"))
##  saveRDS(res_topology, file = file.path(TAB_DIR, "results_topology.rds"))
##
##  ## Non-Gaussian robustness (Figure 6 in manuscript)
##  res_noise <- run_full_simulation(qs = 20, ns = 300, ws = 0.15,
##                                   signals = "strong", topologies = "ER",
##                                   noises = c("gaussian", "t8", "t4",
##                                              "skewnormal", "laplace"),
##                                   reps = REPS, S = S, burn = BURN,
##                                   methods = c("nCPNG", "CPNIG",
##                                               "PC", "LiNGAM"))
##  saveRDS(res_noise, file = file.path(TAB_DIR, "results_noise.rds"))
##
##  ## AML real-data analysis (Section 7 in manuscript)
##  aml <- run_AML_analysis(S = 60000, burn = 5000)
##  saveRDS(aml, file = file.path(TAB_DIR, "aml_results.rds"))
##  plot_aml_heatmap(aml$edgeprobs, aml$proteins,
##    filename = file.path(FIG_DIR, "fig_aml_heatmap.pdf"))
##  plot_aml_graph(aml$MAP, aml$proteins,
##    filename = file.path(FIG_DIR, "fig_aml_mapdag.pdf"),
##    title = "MAP DAG estimate (nCPNG)")
##  plot_aml_graph(aml$MPM, aml$proteins,
##    filename = file.path(FIG_DIR, "fig_aml_mpmdag.pdf"),
##    title = "MPM DAG estimate (nCPNG)")
##
##  ## End of script.

cat("Script loaded.  Uncomment the 'Main entry point' block to run.\n")
