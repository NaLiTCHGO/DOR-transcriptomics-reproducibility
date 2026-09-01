options(save.defaults = list(ascii = FALSE, compress = TRUE, safe = TRUE))
root <- Sys.getenv("DOR_PROJECT_ROOT", unset = ".")
run <- file.path(root, "05_analysis_steps/M02_COHORT_QC/runs/20260813_M02_B1_GSE274832")
out <- file.path(run, "results/expression_qc")

pca <- read.csv(file.path(out, "GSE274832_PCA_COORDINATES.csv"), check.names = FALSE)
cols <- c(DOR = "#D55E00", NOR = "#0072B2")
png(file.path(out, "GSE274832_PCA.png"), width = 1800, height = 1500, res = 240)
par(mar = c(5, 5, 3, 1))
plot(pca$PC1, pca$PC2, pch = 21, bg = cols[pca$phenotype], col = "white", cex = 2.2,
     xlab = sprintf("PC1 (%.1f%%)", pca$PC1_variance_percent[1]),
     ylab = sprintf("PC2 (%.1f%%)", pca$PC2_variance_percent[1]),
     main = "GSE274832: gene-level expression PCA")
text(pca$PC1, pca$PC2, labels = pca$srr, pos = 3, cex = 0.75)
legend("topright", legend = names(cols), pt.bg = cols, pch = 21, bty = "n")
dev.off()

corr <- as.matrix(read.csv(file.path(out, "GSE274832_SAMPLE_CORRELATION.csv"), row.names = 1, check.names = FALSE))
png(file.path(out, "GSE274832_SAMPLE_CORRELATION.png"), width = 1800, height = 1600, res = 240)
par(mar = c(8, 8, 3, 2))
image(seq_len(nrow(corr)), seq_len(ncol(corr)), t(corr[nrow(corr):1, ]),
      col = colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100),
      zlim = c(min(corr), 1), axes = FALSE, main = "GSE274832 sample correlation")
axis(1, at = seq_len(ncol(corr)), labels = colnames(corr), las = 2, cex.axis = 0.75)
axis(2, at = seq_len(nrow(corr)), labels = rev(rownames(corr)), las = 2, cex.axis = 0.75)
box()
dev.off()

markers <- read.csv(file.path(out, "GSE274832_MARKER_EXPRESSION.csv"), check.names = FALSE)
marker_order <- unique(markers$gene_name)
sample_order <- unique(markers$srr)
mat <- xtabs(log2_TPM_plus_1 ~ gene_name + srr, data = markers)[marker_order, sample_order, drop = FALSE]
png(file.path(out, "GSE274832_MARKER_EXPRESSION.png"), width = 1800, height = 2100, res = 240)
par(mar = c(8, 8, 3, 2))
image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, ]),
      col = colorRampPalette(c("white", "#FEE08B", "#D73027"))(100), axes = FALSE,
      main = "Granulosa and retina marker sanity check")
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.75)
axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.75)
box()
dev.off()

cat("PASS: PCA, correlation, and marker QC plots generated\n")
