# Bulk RNA-seq Differential Expression (DESeq2) - R implementation
#
# Which genes respond to SARS-CoV-2 infection in human airway epithelial cells?
# Data: Blanco-Melo et al. 2020, GEO GSE147507 (NHBE, A549, Calu3; mock vs SARS-CoV-2).
# Companion Python version (PyDESeq2) is in rnaseq_de.ipynb - results match.
#
# Run in RStudio: Session -> Set Working Directory -> To Source File Location -> Source.
# First run installs DESeq2, GEOquery, apeglm (Bioconductor) + tidyverse - several minutes.

## ---- Setup ----
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "GEOquery", "apeglm"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(DESeq2); library(GEOquery); library(tidyverse)
dir.create("data_R", showWarnings = FALSE); dir.create("results_R", showWarnings = FALSE)

## ---- Fetch raw counts from GEO ----
if (!file.exists("data_R/GSE147507_RawReadCounts_Human.tsv.gz"))
  getGEOSuppFiles("GSE147507", baseDir = "data_R", makeDirectory = FALSE)
counts <- read.delim(gzfile("data_R/GSE147507_RawReadCounts_Human.tsv.gz"),
                     row.names = 1, check.names = FALSE)

## ---- Metadata + subset to 3 cell lines, mock vs SARS-CoV-2 ----
keep_cols <- grepl("Mock|SARS-CoV-2", colnames(counts)) &
             grepl("NHBE|A549|Calu3",  colnames(counts)) &
             !grepl("ACE2", colnames(counts))
counts <- as.matrix(counts[, keep_cols]); storage.mode(counts) <- "integer"
meta <- data.frame(
  row.names = colnames(counts),
  cell_line = factor(sub("Series[0-9]+_([^_]+)_.*", "\\1", colnames(counts))),
  infection = factor(ifelse(grepl("SARS-CoV-2", colnames(counts)), "SARS_CoV2", "Mock"),
                     levels = c("Mock", "SARS_CoV2"))
)

## ---- DESeq2 with the correct design: control for cell line, test infection ----
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ cell_line + infection)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds)
n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
cat("significant genes (padj < 0.05):", n_sig, "\n")

## ---- The design-formula lesson: why we control for cell line ----
dds_wrong <- DESeq(DESeqDataSetFromMatrix(counts, meta, design = ~ infection)[rowSums(counts) >= 10, ])
cat("~ infection only        -> sig genes:", sum(results(dds_wrong)$padj < 0.05, na.rm = TRUE), "\n")
cat("~ cell_line + infection -> sig genes:", n_sig, "\n")
# Controlling for cell line (which dominates the variance) multiplies detected DE genes.

## ---- Shrink fold-changes (apeglm) for a clean volcano + ranking ----
resLFC <- lfcShrink(dds, coef = "infection_SARS_CoV2_vs_Mock", type = "apeglm")

## ---- Figures ----
# PCA: does variation track infection or cell line?
vsd <- vst(dds, blind = FALSE)
print(plotPCA(vsd, intgroup = c("infection", "cell_line")))
ggsave("results_R/pca.png", width = 7, height = 5, dpi = 150)

# Volcano (shrunken LFC) - hero figure
vdf <- as.data.frame(resLFC)
vdf$significant <- !is.na(vdf$padj) & vdf$padj < 0.05
ggplot(vdf, aes(log2FoldChange, -log10(padj), colour = significant)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c("grey70", "red")) +
  labs(title = "SARS-CoV-2 vs mock - airway epithelium",
       x = "log2 fold change (shrunken)", y = "-log10 adjusted p")
ggsave("results_R/volcano.png", width = 7, height = 5, dpi = 150)

## ---- Save results + top genes ----
res_ordered <- as.data.frame(res)[order(results(dds)$padj), ]
write.csv(res_ordered, "results_R/deseq2_results.csv")
cat("\nTop up-regulated genes on infection:\n")
print(head(res_ordered[!is.na(res_ordered$padj) & res_ordered$padj < 0.05 &
                       res_ordered$log2FoldChange > 0, ], 15))

## ---- Functional enrichment: GO/KEGG over-representation + GSEA ----
# Turn the DE gene list into pathways. ORA = hypergeometric test of the significant UP genes
# against GO/KEGG (background = all tested genes); GSEA = rank ALL genes and find coordinately
# shifted pathways (NES > 0 up, < 0 down). clusterProfiler needs Entrez IDs, so map symbols first.
for (p in c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
de <- as.data.frame(res); de <- de[!is.na(de$padj) & !is.na(de$log2FoldChange), ]
ids <- bitr(rownames(de), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de$SYMBOL <- rownames(de); de <- merge(de, ids, by = "SYMBOL")

sig_up <- de$ENTREZID[de$padj < 0.05 & de$log2FoldChange > 1]                 # ORA input
ego <- enrichGO(sig_up, universe = de$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
write.csv(as.data.frame(ego), "results_R/enrichGO_up.csv", row.names = FALSE)
if (nrow(as.data.frame(ego)) > 0)
  ggsave("results_R/enrichGO_up_dotplot.png", dotplot(ego, showCategory = 15), width = 8, height = 6, dpi = 150)

ekegg <- enrichKEGG(sig_up, organism = "hsa", pvalueCutoff = 0.05)
write.csv(as.data.frame(ekegg), "results_R/enrichKEGG_up.csv", row.names = FALSE)

geneList <- sort(setNames(sign(de$log2FoldChange) * -log10(de$pvalue + 1e-300), de$ENTREZID),
                 decreasing = TRUE)                                           # GSEA: rank all genes
gse <- gseKEGG(geneList, organism = "hsa", pvalueCutoff = 0.05)
write.csv(as.data.frame(gse), "results_R/gseKEGG.csv", row.names = FALSE)
if (nrow(as.data.frame(gse)) > 0)
  ggsave("results_R/gsea_top.png", gseaplot2(gse, geneSetID = 1), width = 8, height = 6, dpi = 150)
cat("\nTop GSEA KEGG (by NES):\n"); print(head(as.data.frame(gse)[, c("Description", "NES", "p.adjust")], 10))

## ---- Interpretation ----
# PCA is dominated by cell line (PC1 70% + PC2 24% = ~94% of variance); infection is a
# subtle within-cell-line shift. That is exactly why the design formula matters: with
# ~ infection alone the cell-line variance swamps the signal (873 sig genes), while
# ~ cell_line + infection removes it and reveals the infection effect (3976 sig genes).
# Top up-regulated genes are inflammatory cytokines (IL36G, IL1A) and antiviral ISGs
# (MX1) - the imbalanced cytokine/interferon response Blanco-Melo et al. reported.
#
# Enrichment collapses that gene list into one coherent program: ORA's top GO term is
# cytokine-mediated signaling (adj p ~1e-26), and GSEA ranks cytokine-cytokine receptor
# interaction, TNF, NF-kappa B, NOD-like receptor and JAK-STAT signalling at the top
# (NES ~ +2.4-2.6, FDR q ~ 0), while oxidative phosphorylation is the strongest DOWN
# pathway (NES ~ -2.5). Biologically: infection drives a coordinated innate-immune /
# inflammatory response and suppresses host energy metabolism - a pathway-level result,
# not a random gene list. Caveats: annotation bias; ORA depends on threshold + universe;
# pathways overlap; enrichment is correlational.
