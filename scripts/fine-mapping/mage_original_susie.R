# Model and credible-set export follow MAGE run_susie.R (mccoy-lab/MAGE, f929fb5).
suppressPackageStartupMessages(library(susieR))
stopifnot(as.character(packageVersion("susieR")) == "0.12.16")
set.seed(1)
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
samples <- read.delim(file.path(input, "samples.tsv"), check.names = FALSE)
variants <- read.delim(file.path(input, "z.tsv"), stringsAsFactors = FALSE)
n <- nrow(samples)
p <- nrow(variants)
X <- matrix(readBin(file.path(input, "dosage.bin"), integer(), n * p,
                    size = 1, signed = FALSE), nrow = n, ncol = p)
Z <- cbind(1, as.matrix(samples[, -(1:2)]))
y <- samples$expression
stopifnot(n == 731, all(X %in% 0:2), all(is.finite(variants$z)))
A <- crossprod(Z)
X <- X - Z %*% solve(A, crossprod(Z, X))
y <- y - Z %*% solve(A, crossprod(Z, y))
# BLAS cross-product of standardized columns is the sample correlation matrix.
X_scaled <- scale(X)
marginal_r <- as.vector(crossprod(X_scaled, scale(y))) / (n - 1)
z_from_data <- marginal_r * sqrt((n - qr(Z)$rank - 1) / (1 - marginal_r^2))
LD <- crossprod(X_scaled) / (n - 1)
rm(X_scaled, X)
fit <- susie_rss(z = variants$z, n = n, R = LD, var_y = var(as.vector(y)),
                 L = 10, coverage = 0.95, min_abs_corr = 0.5,
                 estimate_residual_variance = TRUE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(snp_id = variants$snp_id, pip = fit$pip),
            file.path(output, "variants.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# The published exporter assigns one label per SNP, overwriting overlaps in set order.
component <- rep(NA_character_, p)
for (label in names(fit$sets$cs)) component[fit$sets$cs[[label]]] <- label
sets <- na.omit(data.frame(component = component, snp_id = variants$snp_id))
write.table(sets, file.path(output, "credible_sets.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(n = n, variants = p, iterations = fit$niter, converged = fit$converged,
                       nominal_z_correlation = cor(z_from_data, variants$z)),
            file.path(output, "fit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
