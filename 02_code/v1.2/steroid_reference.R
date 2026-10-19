## v1.2 plan section 2: in-vitro glucocorticoid (GC) reference set. No cardiac data are read in this script.
## Per dataset: limma (voom for counts; limma-trend for normalised RNA-seq; lmFit for arrays; one-sample for two-colour log-ratios),
## donor-paired where the design allows. QC: FKBP5 and TSC22D3 logFC>1 & BH-FDR<0.05. Meta: DerSimonian-Laird random effects
## on logFC (SE = sqrt(s2.post)*stdev.unscaled), genes measured in >=2 QC-passing datasets. Classification per plan.
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(readxl)})
setwd("v12/steroid/data"); set.seed(1)
res <- list(); info <- list()
cc <- function(E) E[apply(E, 1, function(r) all(is.finite(r))), , drop = FALSE]  # drop rows with NA/Inf (code fix 2026-10-03)
fit_out <- function(fit, coef, name, design_note) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  d <- data.table(gene = rownames(tt), logFC = tt$logFC, se = as.numeric(se), t = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val, AveExpr = tt$AveExpr)
  info[[name]] <<- design_note; res[[name]] <<- d; invisible(d)
}
collapse_max <- function(E, sym) {   # keep, per symbol, the row with the highest mean expression
  ok <- !is.na(sym) & sym != "" & !grepl("///", sym); E <- E[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(sym, -rowMeans(E)); E <- E[o, , drop = FALSE]; sym <- sym[o]; k <- !duplicated(sym); E <- E[k, , drop = FALSE]; rownames(E) <- sym[k]; E
}
read_sm <- function(f) {  # GEO series matrix -> list(expr, sample titles, characteristics)
  L <- readLines(gzfile(f)); b <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
  X <- fread(text = paste(L[(b + 1):(e - 1)], collapse = "\n"), header = TRUE); m <- as.matrix(X[, -1]); rownames(m) <- X[[1]]
  tit <- gsub('"', "", strsplit(grep("^!Sample_title", L, value = TRUE), "\t")[[1]][-1])
  ch <- lapply(grep("^!Sample_characteristics_ch1", L, value = TRUE), function(z) gsub('"', "", strsplit(z, "\t")[[1]][-1]))
  list(E = m, title = tit, ch = ch)
}
## 1 GSE150049: DESeq-normalised counts; pairing not documented -> unpaired; limma-trend on log2(x+1)
x <- as.data.frame(read_excel("GSE150049_counts.xlsx")); x <- x[!duplicated(x$GeneName) & !is.na(x$GeneName), ]
m <- as.matrix(x[, -1]); rownames(m) <- x$GeneName; grp <- factor(ifelse(as.integer(sub("Mono", "", colnames(m))) >= 51, "GC", "ctrl"), levels = c("ctrl", "GC"))
lm_ <- log2(m + 1); keep <- rowSums(m >= 1) >= min(table(grp)); lm_ <- lm_[keep, ]
fit <- eBayes(lmFit(cc(lm_), model.matrix(~grp)), trend = TRUE); fit_out(fit, "grpGC", "GSE150049", sprintf("monocytes dex vs untreated, unpaired, %d v %d", sum(grp == "GC"), sum(grp == "ctrl")))
## 2 GSE85333: MDM prednisolone vs vehicle CTRL, unstimulated 24 h, 6 donors paired (values linear -> log2)
s <- read_sm("GSE85333_series_matrix.txt.gz"); chs <- sapply(s$ch, function(z) sub("^[^:]+: ", "", z))
colnames(chs) <- sapply(s$ch, function(z) sub(":.*", "", z[1])); stim <- chs[, grep("inflammatory", colnames(chs))]; drug <- chs[, grep("drug", colnames(chs))]; donor <- chs[, "donor"]
k <- grepl("No inflammation", stim) & drug %in% c("2 CTRL", "7 Prednisolone"); E <- log2(pmax(s$E[, k], 1)); g <- factor(ifelse(drug[k] == "7 Prednisolone", "GC", "ctrl"), levels = c("ctrl", "GC")); dn <- factor(donor[k])
E <- cc(E); E <- E[rowMeans(E) > log2(10), ]
fit <- eBayes(lmFit(E, model.matrix(~dn + g)), trend = TRUE); fit_out(fit, "gGC", "GSE85333", sprintf("MDM prednisolone vs vehicle (unstimulated 24 h), paired, %d donors", nlevels(dn)))
## 3 GSE35683: cord-blood monocytes LPS+dex vs LPS+PBS, 5 donors paired by index (GPL570, log2)
s <- read_sm("GSE35683_series_matrix.txt.gz"); tt <- s$title; k <- grepl("^Mono-(DEX|LPS)-", tt)
gpl <- fread("GPL570.txt", skip = "ID\tGB_ACC", sep = "\t", quote = "", fill = TRUE, select = c("ID", "Gene Symbol"))
sym <- gpl$`Gene Symbol`[match(toupper(rownames(s$E)), toupper(gpl$ID))]
E <- collapse_max(s$E[, k], sym); g <- factor(ifelse(grepl("DEX", tt[k]), "GC", "ctrl"), levels = c("ctrl", "GC")); dn <- factor(sub("^Mono-(DEX|LPS)-(\\d)-.*", "\\2", tt[k]))
E <- cc(E); E <- E[rowMeans(E) > quantile(rowMeans(E), 0.25), ]
fit <- eBayes(lmFit(E, model.matrix(~dn + g)), trend = TRUE); fit_out(fit, "gGC", "GSE35683", "cord-blood monocytes LPS+dex vs LPS+PBS (stimulated design), paired, 5 donors")
## 4 GSE112101: monocytes methylprednisolone vs vehicle, 6 h, paired BR1-4; edgeR TMM + voom
x <- fread("GSE112101_RawReadCounts.txt.gz"); m <- as.matrix(x[, -1]); rownames(m) <- x$gene; m <- m[!duplicated(rownames(m)), ]
k <- grepl("^Mono_6H_(MP|VH)_BR", colnames(m)); m <- m[, k]; g <- factor(ifelse(grepl("_MP_", colnames(m)), "GC", "ctrl"), levels = c("ctrl", "GC")); dn <- factor(sub(".*_BR", "", colnames(m)))
y <- DGEList(m); des <- model.matrix(~dn + g); y <- y[filterByExpr(y, des), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
fit <- eBayes(lmFit(voom(y, des), des)); fit_out(fit, "gGC", "GSE112101", "monocytes methylprednisolone vs vehicle 6 h, paired, 4 donors")
## 5 E-MEXP-2598: monocytes dex 36 h vs no dex; 4 direct two-colour arrays, Cy5 = dex on all arrays -> one-sample on log2 ratio
arr <- fread("EMEXP2598_dex0_arrays.csv"); a <- unique(arr[, .(f = `Array Data File`, dexlab = Label[grepl("^Dexa", `Source Name`)]), by = `Hybridization Name`])
stopifnot(all(a$dexlab == "Cy5"))
M <- NULL; A <- NULL
for (f in a$f) {
  z <- fread(f, sep = "\t", select = c("Feature Extraction Software:ControlType", "Feature Extraction Software:GeneName", "Feature Extraction Software:LogRatio",
                                         "Feature Extraction Software:gProcessedSignal", "Feature Extraction Software:rProcessedSignal",
                                         "Feature Extraction Software:gIsWellAboveBG", "Feature Extraction Software:rIsWellAboveBG"))
  setnames(z, c("ctl", "gene", "lr", "g", "r", "gw", "rw")); z <- z[ctl == 0 & (gw == 1 | rw == 1) & !is.na(lr)]
  z <- z[, .(lr = mean(lr * log2(10)), a = mean(0.5 * log2(pmax(g, 1)) + 0.5 * log2(pmax(r, 1)))), by = gene]
  M <- if (is.null(M)) z[, .(gene, lr)] else merge(M, z[, .(gene, lr)], by = "gene", suffixes = c("", f))
  A <- if (is.null(A)) z[, .(gene, a)] else merge(A, z[, .(gene, a)], by = "gene", suffixes = c("", f))
}
Mm <- as.matrix(M[, -1]); rownames(Mm) <- M$gene; Am <- rowMeans(as.matrix(A[, -1])); keepA <- Am > quantile(Am, 0.25, na.rm = TRUE) & !is.na(Am); Mm <- cc(Mm[keepA, , drop = FALSE])
fit <- eBayes(lmFit(Mm, matrix(1, ncol(Mm), 1, dimnames = list(NULL, "GC"))), trend = FALSE); fit_out(fit, "GC", "E-MEXP-2598", "monocytes dex 36 h vs no dex, two-colour direct, 4 donors")
## 6 E-MTAB-7572: MDM dex 24 h vs no stimulus (M0), paired d54/d55/d57; entrez counts -> symbols; voom
x <- fread("EMTAB7572_entrez_counts.csv"); gi <- fread("Homo_sapiens.gene_info.gz", select = c("GeneID", "Symbol"))
m <- as.matrix(x[, -1]); sym <- gi$Symbol[match(x[[1]], gi$GeneID)]; rownames(m) <- sym; m <- m[!is.na(sym) & !duplicated(sym), ]
k <- grepl("d5[457](D|M0)_", colnames(m)); m <- m[, k]; g <- factor(ifelse(grepl("D_S", colnames(m)), "GC", "ctrl"), levels = c("ctrl", "GC")); dn <- factor(sub(".*(d5[457]).*", "\\1", colnames(m)))
stopifnot(sum(g == "GC") == 3, sum(g == "ctrl") == 3)
y <- DGEList(m); des <- model.matrix(~dn + g); y <- y[filterByExpr(y, des), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
fit <- eBayes(lmFit(voom(y, des), des)); fit_out(fit, "gGC", "E-MTAB-7572", "MDM dex 24 h vs unstimulated, paired, 3 donors")
## ---- QC ----
qc <- rbindlist(lapply(names(res), function(n) { d <- res[[n]]; q <- d[gene %in% c("FKBP5", "TSC22D3")]
  data.table(dataset = n, design = info[[n]], n_genes = nrow(d), FKBP5_logFC = q[gene == "FKBP5", logFC], FKBP5_FDR = q[gene == "FKBP5", FDR],
             TSC22D3_logFC = q[gene == "TSC22D3", logFC], TSC22D3_FDR = q[gene == "TSC22D3", FDR]) }), fill = TRUE)
qc[, pass := !is.na(FKBP5_logFC) & !is.na(TSC22D3_logFC) & FKBP5_logFC > 1 & FKBP5_FDR < 0.05 & TSC22D3_logFC > 1 & TSC22D3_FDR < 0.05]
fwrite(qc, "../dataset_QC.csv"); print(qc)
keepds <- qc[pass == TRUE, dataset]; stopifnot(length(keepds) >= 2)
long <- rbindlist(lapply(names(res), function(n) res[[n]][, dataset := n])); fwrite(long, "../per_dataset_log2fc.csv.gz")
## ---- DL random-effects meta (QC-passing datasets only) ----
L <- long[dataset %in% keepds & is.finite(se) & se > 0]
mt <- L[, { w <- 1 / se^2; k <- .N; bf <- sum(w * logFC) / sum(w); Q <- sum(w * (logFC - bf)^2)
  tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0; ws <- 1 / (se^2 + tau2); mu <- sum(ws * logFC) / sum(ws); sem <- sqrt(1 / sum(ws))
  .(k = k, meta_log2FC = mu, se = sem, ci_lo = mu - 1.96 * sem, ci_hi = mu + 1.96 * sem, P = 2 * pnorm(-abs(mu / sem)), tau2 = tau2,
    n_pos = sum(logFC > 0), n_neg = sum(logFC < 0)) }, by = gene][k >= 2]
mt[, same_sign_frac := pmax(n_pos, n_neg) / k]
mt[, class := fifelse(meta_log2FC >= 0.5 & P < 0.01 & n_pos / k >= 2/3, "GC_induced",
               fifelse(meta_log2FC <= -0.5 & P < 0.01 & n_neg / k >= 2/3, "GC_repressed",
               fifelse(ci_lo > -0.5 & ci_hi < 0.5, "GC_insensitive", "ambiguous")))]
fwrite(mt, "../meta_all_genes.csv.gz")
mods <- fread("../../../modules_v1.1.tsv")[module %in% c("RES_EFF_core", "RES_EFF_lit", "MHCII_core", "ISG_IFN", "MAC_PAN_ID"), .(module, gene)]
mods <- rbind(mods, data.table(module = "GC_response", gene = c("FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A")))
cls <- merge(mods, mt, by = "gene", all.x = TRUE); cls[is.na(class), class := "not_measured_in_>=2"]
fwrite(cls[order(module, class, gene)], "../gc_gene_classification.csv")
print(cls[, .N, by = .(module, class)][order(module, class)])
print(cls[module == "RES_EFF_core", .(gene, k, meta_log2FC = round(meta_log2FC, 2), ci_lo = round(ci_lo, 2), ci_hi = round(ci_hi, 2), P = signif(P, 2), same_sign_frac = round(same_sign_frac, 2), class)])
