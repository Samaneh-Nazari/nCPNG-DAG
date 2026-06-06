## ============================================================================
##  Real-data analysis (AML protein expression) for the JCGS revision of
##  "Bayesian Structure Learning of Gaussian DAGs via Non-conjugate
##   Normal-Gamma Priors with Adaptive MCMC"
##
##  This script loads the simulation/methodological functions from
##  nCPNG_simulation.R and applies them to the leukaemia protein
##  expression data of Kornblau et al., as redistributed in the
##  BCDAG R package.
##
##  Outputs:
##    - figs/fig_mcmc_diagnostics.pdf  (4-chain trace and R-hat)
##    - figs/fig_aml_trace.pdf         (graph-size trace + running mean)
##    - figs/fig_aml_heatmap.pdf       (posterior edge probabilities)
##    - figs/fig_aml_mapdag.pdf        (MAP DAG)
##    - figs/fig_aml_mpmdag.pdf        (MPM DAG)
##    - tables/aml_posterior_summary.csv
## ============================================================================

source("nCPNG_simulation.R")

## ----- Load data ------------------------------------------------------------
load_AML <- function() {
  ## Preferred: from the BCDAG package
  if (requireNamespace("BCDAG", quietly = TRUE)) {
    return(as.matrix(BCDAG::leukemia))
  }
  ## Fallback: local copy
  candidates <- c("leukemia.rda", "data/leukemia.rda")
  for (f in candidates) {
    if (file.exists(f)) { load(f); return(as.matrix(leukemia)) }
  }
  stop("Cannot locate the AML dataset.  Install BCDAG or supply leukemia.rda.")
}

X <- load_AML()
n <- nrow(X); q <- ncol(X)
proteins <- colnames(X)
if (is.null(proteins)) proteins <- paste0("P", seq_len(q))
cat(sprintf("Loaded AML data: n=%d patients, q=%d proteins\n", n, q))

## ----- Hyperparameter settings (manuscript Section 7) -----------------------
ALPHA <- q
U     <- diag(1 / n, q)
W     <- 0.5
S     <- 60000
BURN  <- 5000

## ----- 4-chain diagnostics --------------------------------------------------
cat("Running 4-chain MCMC for Gelman-Rubin diagnostics...\n")
diag_out <- multi_chain_diagnostics(X, n_chains = 4, S = S, burn = BURN,
                                     a = ALPHA, U = U, w = W)
cat(sprintf("Gelman-Rubin R-hat (graph size): %.4f\n", diag_out$R.hat))
cat(sprintf("Total ESS across chains: %.0f\n", diag_out$ESS))

## ----- Main posterior run ---------------------------------------------------
cat("Running main posterior chain...\n")
set.seed(1)
out <- learn_DAG_nCPNG(S = S, burn = BURN, data = X,
                       a = ALPHA, U = U, w = W, fast = TRUE)

## ----- Posterior summaries --------------------------------------------------
edgeprobs <- apply(out$Graphs, c(1, 2), mean)
dimnames(edgeprobs) <- list(proteins, proteins)
MPM <- get_MPMdag(out); dimnames(MPM) <- list(proteins, proteins)
MAP <- get_MAPdag(out); dimnames(MAP) <- list(proteins, proteins)
graph.size <- out$graph.size

post_mean_size <- mean(graph.size)
post_sd_size   <- sd(graph.size)
n_MPM <- sum(MPM)
n_MAP <- sum(MAP)

cat(sprintf("Posterior mean graph size: %.2f  (SD %.2f)\n",
            post_mean_size, post_sd_size))
cat(sprintf("MPM edges: %d   |   MAP edges: %d\n", n_MPM, n_MAP))

## ----- Write CSV summary ----------------------------------------------------
dir.create("tables", showWarnings = FALSE)
write.csv(round(edgeprobs, 4),
          file = "tables/aml_edge_inclusion_probabilities.csv")
write.csv(MPM, file = "tables/aml_MPM_adjacency.csv")
write.csv(MAP, file = "tables/aml_MAP_adjacency.csv")

summary_df <- data.frame(
  quantity = c("Posterior mean graph size", "Posterior SD graph size",
               "Number of edges in MPM",  "Number of edges in MAP",
               "Gelman-Rubin R-hat", "Total ESS"),
  value    = c(post_mean_size, post_sd_size, n_MPM, n_MAP,
               diag_out$R.hat, diag_out$ESS))
write.csv(summary_df, file = "tables/aml_posterior_summary.csv",
          row.names = FALSE)

## ----- Figures --------------------------------------------------------------
dir.create("figs", showWarnings = FALSE)
plot_diagnostics(diag_out, filename = "figs/fig_mcmc_diagnostics.pdf")

##  Trace + running mean
pdf("figs/fig_aml_trace.pdf", width = 11, height = 3.4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), family = "serif")
plot(seq_along(graph.size), graph.size, type = "l",
     col = "#1F4E79", lwd = 0.4,
     xlab = "Iteration", ylab = "Graph size",
     main = "(a) MCMC trace of graph size (nCPNG)")
running_mean <- cumsum(graph.size) / seq_along(graph.size)
plot(seq_along(graph.size), running_mean, type = "l",
     col = "#1F4E79", lwd = 1.2,
     xlab = "Iteration", ylab = "Running posterior mean",
     main = "(b) Running posterior mean")
dev.off()

plot_aml_heatmap(edgeprobs, proteins,
                 filename = "figs/fig_aml_heatmap.pdf")
plot_aml_graph(MAP, proteins,
               filename = "figs/fig_aml_mapdag.pdf",
               title = "MAP DAG estimate (nCPNG)")
plot_aml_graph(MPM, proteins,
               filename = "figs/fig_aml_mpmdag.pdf",
               title = "MPM DAG estimate (nCPNG)")

cat("\nAll AML outputs written to figs/ and tables/.\n")
