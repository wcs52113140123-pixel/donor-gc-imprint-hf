# gc_tools.R — R port of gc_tools.py (same definitions; see header of gc_tools.py and the v1.3 report)
GC_SPECIFIC3 <- c("FKBP5", "TSC22D3", "ZBTB16")
GC_RESPONSE8 <- c("FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A")
gc_load_table <- function(path = "gc_gene_classification_all.csv") read.csv(path)

# Z: samples x genes matrix of reference-standardised expression, colnames = gene symbols
gc_score <- function(Z, ref, which = c("GC_specific3", "GC_response")) {
  gl <- if (match.arg(which) == "GC_specific3") GC_SPECIFIC3 else GC_RESPONSE8
  g <- intersect(gl, colnames(Z)); if (length(g) < 2) stop("fewer than 2 score genes measured")
  s <- rowMeans(Z[, g, drop = FALSE], na.rm = TRUE); (s - mean(s[ref])) / sd(s[ref])
}
gc_flag <- function(df, gene_col = "gene", table = gc_load_table()) {
  t <- table[, c("gene", "class", "meta_log2FC")]; names(t) <- c(gene_col, "GC_class", "GC_meta_log2FC")
  out <- merge(df, t, by = gene_col, all.x = TRUE, sort = FALSE); out$GC_class[is.na(out$GC_class)] <- "not_in_reference"; out
}
# design: model matrix with intercept in column 1 and the contrast of interest in column `coef` (default 2)
gc_collinearity <- function(design, score, coef = 2) {
  D0 <- as.matrix(design); D1 <- cbind(D0, score); cs <- D0[, coef] == 1
  vif <- MASS::ginv(crossprod(D1))[coef, coef] / MASS::ginv(crossprod(D0))[coef, coef]
  auc <- unname(wilcox.test(score[!cs], score[cs], exact = FALSE)$statistic) / (sum(!cs) * sum(cs))
  y <- D0[, coef]; r2 <- 1 - sum(resid(lm(y ~ D1[, -coef] - 1))^2) / sum((y - mean(y))^2)
  list(AUC = auc, R2_case = r2, VIF = vif, identifiable = (auc > 0.1 & auc < 0.9 & vif < 5))
}
# Adds the GC score as a covariate; with limma when available (recommended for counts -> voom), else OLS.
gc_adjusted_de <- function(Y, design, score, coef = 2, check = TRUE) {
  if (check) { cc <- gc_collinearity(design, score, coef); if (!cc$identifiable) warning("GC score nearly aliased with the contrast; adjusted DE not interpretable") }
  D <- cbind(as.matrix(design), GC = score)
  fit <- limma::eBayes(limma::lmFit(Y, D)); limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
}
