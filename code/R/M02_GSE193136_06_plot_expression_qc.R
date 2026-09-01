options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- if (length(script_arg)) normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "../.."), winslash = "/") else getwd()
root <- Sys.getenv("DOR_PROJECT_ROOT", unset = default_root)
run <- file.path(root, "05_analysis_steps/M02_COHORT_QC/runs/20260814_M02_B3_GSE193136")
out <- file.path(run, "results/expression_qc")

pca <- read.csv(file.path(out, "GSE193136_PCA_COORDINATES.csv"), check.names = FALSE)
cols <- c(DOR = "#D55E00", NOR = "#0072B2")
png(file.path(out, "GSE193136_PCA.png"), width = 2200, height = 1600, res = 240)
par(mar = c(5, 5, 3, 7), xpd = NA)
xpad <- diff(range(pca$PC1)) * 0.12
ypad <- diff(range(pca$PC2)) * 0.12
plot(pca$PC1, pca$PC2, pch = 21, bg = cols[pca$phenotype], col = "white", cex = 2.0,
     xlim = range(pca$PC1) + c(-xpad, xpad), ylim = range(pca$PC2) + c(-ypad, ypad),
     xlab = sprintf("PC1 (%.1f%%)", pca$PC1_variance_percent[1]),
     ylab = sprintf("PC2 (%.1f%%)", pca$PC2_variance_percent[1]),
     main = "GSE193136: gene-level expression PCA")
label_pos <- c(3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 4, 1)
text(pca$PC1, pca$PC2, labels = paste0(pca$srr, " (", pca$age_years, ")"), pos = label_pos, cex = 0.66)
legend("topright", inset = c(-0.16, 0), legend = names(cols), pt.bg = cols, pch = 21, bty = "n")
dev.off()

corr <- as.matrix(read.csv(file.path(out, "GSE193136_SAMPLE_CORRELATION.csv"), row.names = 1, check.names = FALSE))
png(file.path(out, "GSE193136_SAMPLE_CORRELATION.png"), width = 2100, height = 1900, res = 240)
par(mar = c(9, 9, 3, 2))
image(seq_len(nrow(corr)), seq_len(ncol(corr)), t(corr[nrow(corr):1, ]),
      col = colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100),
      zlim = c(min(corr), 1), axes = FALSE, xlab = "", ylab = "",
      main = "GSE193136 sample correlation")
axis(1, at = seq_len(ncol(corr)), labels = colnames(corr), las = 2, cex.axis = 0.65)
axis(2, at = seq_len(nrow(corr)), labels = rev(rownames(corr)), las = 2, cex.axis = 0.65)
box(); dev.off()

markers <- read.csv(file.path(out, "GSE193136_MARKER_EXPRESSION.csv"), check.names = FALSE)
marker_order <- unique(markers$gene_name)
sample_order <- unique(markers$srr)
mat <- xtabs(log2_TPM_plus_1 ~ gene_name + srr, data = markers)[marker_order, sample_order, drop = FALSE]
png(file.path(out, "GSE193136_MARKER_EXPRESSION.png"), width = 2200, height = 2400, res = 240)
par(mar = c(9, 9, 3, 2))
image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, ]),
      col = colorRampPalette(c("white", "#FEE08B", "#D73027"))(100), axes = FALSE,
      xlab = "", ylab = "",
      main = "Granulosa and off-target marker sanity check")
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.62)
axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.66)
box(); dev.off()

writeLines("PASS\tPCA, correlation and marker QC plots generated", file.path(run, "STEP06_QC_PLOTS.PASS.txt"))
cat("PASS: GSE193136 PCA, correlation and marker QC plots generated\n")
