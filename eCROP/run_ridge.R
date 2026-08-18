#!/usr/bin/env Rscript
# Ridge regression on scImpute data (standalone)
suppressPackageStartupMessages({
  library(Seurat); library(glmnet); library(Matrix)
})

# 请根据实际环境修改以下路径
output_dir <- "/path/to/project/res/eCROP/scimpute"
ridge_dir <- file.path(output_dir, "ridge_output")
dir.create(ridge_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Ridge regression ===\n")

# Load full data
obj <- readRDS(file.path(output_dir, "seurat_with_imputed.rds"))
ep_all <- as.matrix(GetAssayData(obj, assay = "scImpute", layer = "data"))

# Assign gene from STC_t
stc_t <- GetAssayData(obj, assay = "STC_t", layer = "data")
assign_gene <- apply(stc_t, 2, function(x) {
  mx <- which.max(x)
  if (x[mx] > 0) rownames(stc_t)[mx] else "NA"
})

keep <- assign_gene != "NA"
ep <- ep_all[, keep]
pi <- assign_gene[keep]

# Use top 1000 HVGs
ep_var <- apply(ep, 1, var)
top_genes <- names(sort(ep_var, decreasing = TRUE))[1:1000]
ep <- ep[top_genes, ]
cat(sprintf("  Expression: %d genes x %d cells\n", nrow(ep), ncol(ep)))

# Build one-hot design
all_perturbs <- sort(unique(pi))
ko_levels <- all_perturbs[all_perturbs != "STC"]
ko_cell_counts <- table(pi[pi != "STC"])
valid_kos <- names(ko_cell_counts[ko_cell_counts >= 5])
keep2 <- pi %in% c("STC", valid_kos)
ep <- ep[, keep2]
pi <- pi[keep2]
ko_levels <- valid_kos
cat(sprintf("  KOs: %d, Cells: %d\n", length(ko_levels), ncol(ep)))

n_cells <- ncol(ep)
n_features <- length(ko_levels)
X <- Matrix(0, nrow = n_cells, ncol = n_features,
            dimnames = list(colnames(ep), ko_levels), sparse = TRUE)
for (k in ko_levels) X[pi == k, k] <- 1

all_targets <- rownames(ep)
n_targets <- length(all_targets)
coef_matrix <- matrix(0, nrow = n_features, ncol = n_targets,
                      dimnames = list(ko_levels, all_targets))

cat(sprintf("  Training %d Ridge models...\n", n_targets))
for (i in seq_len(n_targets)) {
  y <- as.numeric(ep[i, ])
  fit <- glmnet(X, y, alpha = 0, lambda = 0.1, standardize = FALSE, intercept = TRUE)
  coef_matrix[, i] <- as.numeric(coef(fit)[-1])
  if (i %% 100 == 0) cat(sprintf("    %d/%d\n", i, n_targets))
}

write.csv(coef_matrix, file.path(ridge_dir, "coef_matrix_all.csv"))
cat(sprintf("  Full coefficient matrix: %d x %d\n", nrow(coef_matrix), ncol(coef_matrix)))

# TF-TF subnetwork
tf_coef <- coef_matrix[, intersect(colnames(coef_matrix), rownames(coef_matrix))]
tf_coef <- tf_coef[intersect(rownames(tf_coef), colnames(tf_coef)),
                    intersect(colnames(tf_coef), rownames(tf_coef))]
write.csv(tf_coef, file.path(ridge_dir, "tf_coefficient_matrix.csv"))
cat(sprintf("  TF-TF subnetwork: %d x %d\n", nrow(tf_coef), ncol(tf_coef)))

cat("Ridge regression done.\n")
