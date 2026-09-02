# Bulk RNA-seq Differential Expression — SARS-CoV-2 Host Response

Which genes respond to **SARS-CoV-2 infection** in human airway epithelial cells? A differential-expression
(DE) analysis of Blanco-Melo et al. 2020 (**GEO GSE147507**), implemented in **R (DESeq2)** and **Python
(PyDESeq2)** with matching results. The project doubles as a concrete lesson in why the **design formula** —
not just the code — decides whether a DE result is right.

![Volcano: SARS-CoV-2 vs mock](results_R/volcano.png)

## 1. Research question

Which genes are differentially expressed when human airway epithelial cells are infected with SARS-CoV-2
versus mock — and how much does controlling for the cell-line confounder change the answer?

## 2. Biological background

**RNA-seq** yields a **count matrix** (genes × samples). Counts are non-negative integers whose variance
grows with the mean, so DE tools model them with the **negative binomial** distribution (not a t-test),
after normalising for sequencing depth via **size factors**. For each gene we get a **log₂ fold-change** and
a **p-value**; testing ~20k genes requires **FDR-adjusted p-values** (`padj`, Benjamini–Hochberg), and we
call genes significant at `padj < 0.05`. The data: primary bronchial (NHBE), A549, and Calu3 cells, each
mock- or SARS-CoV-2-infected.

## 3. Dataset

- **Source:** GEO **GSE147507** (Blanco-Melo et al., *Cell* 2020), human raw read counts.
- **Access:** fetched programmatically — `GEOquery::getGEOSuppFiles()` (R) / pandas from the GEO FTP (Python).
- **Samples used:** NHBE, A549, Calu3 — mock vs SARS-CoV-2 (~9 vs 9 across three cell lines).

## 4. Methods

DESeq2 / PyDESeq2 with design **`~ cell_line + infection`** — control for cell line, then test infection.
The pipeline: size factors → dispersion estimation → negative-binomial Wald test → FDR. Fold-changes are
shrunk with **`apeglm`** (`lfcShrink`) for a clean volcano and reliable ranking. Sample names are parsed into
a metadata table (cell line + infection status) — a routine real-world data-wrangling step.

## 5. Analysis pipeline

1. Fetch the raw count matrix from GEO.
2. Parse sample names → metadata; subset to the three cell lines, mock vs SARS-CoV-2.
3. Build the DESeqDataSet with `~ cell_line + infection`; pre-filter low-count genes.
4. Run DESeq2; extract results (log₂FC, `padj`); count significant genes.
5. Compare designs (`~ infection` vs `~ cell_line + infection`) to show the confounder effect.
6. Shrink fold-changes; plot PCA + volcano; save results.

## 6. Results

| metric | value |
|--------|-------|
| significant genes, `~ cell_line + infection` (`padj<0.05`) | **3,976** |
| significant genes, `~ infection` only | 873 |
| effect of controlling for cell line | **4.5× more DE genes** |
| top up-regulated genes | IL36G, IL1A, SOD2, IRAK2, BIRC3, **MX1** |

Full results in [`results_R/deseq2_results.csv`](results_R/deseq2_results.csv). R and Python give identical
top genes and fold-changes (e.g. IL36G log₂FC ≈ 2.37, IL1A ≈ 3.25 in both).

## 7. Figures

`results_R/volcano.png` — shrunken-LFC volcano (hero, above); `pca.png` — sample PCA. Python equivalents in
`results_py/`.

## 8. Limitations

- Three cell lines with few replicates each — modest power; results are cell-culture, not in-vivo, biology.
- Genomic-span/count-level analysis only; no isoform- or allele-level resolution.
- No pathway/enrichment step here (a natural extension); interpretation is based on top individual genes.
- DE detects association with infection status, not causal/mechanistic direction.

## 9. Biological interpretation

The **PCA is dominated by cell line** (PC1 70% + PC2 24% ≈ 94% of variance); the infection effect is a subtle
within-cell-line shift. This is precisely why the **design formula** matters: with `~ infection` alone the
cell-line variance swamps the signal (873 genes), whereas `~ cell_line + infection` removes it and reveals
the infection response (**3,976 genes**). The top up-regulated genes are **inflammatory cytokines** (IL36G,
IL1A) and **antiviral interferon-stimulated genes** (MX1) — the imbalanced cytokine/interferon host response
to SARS-CoV-2 reported by Blanco-Melo et al.

## Language comparison (R vs Python)

Implemented in **DESeq2** (R) and **PyDESeq2** (Python) on the same data and design; the two agree to the
decimal on top genes and fold-changes — a strong cross-toolchain validation. The main idiomatic difference:
**DESeq2 expects counts as genes × samples, PyDESeq2 as samples × genes** (a transpose) — the most common
gotcha when porting a DE analysis between the two.

## 10. Reproduce

```bash
# R (RStudio)
install.packages(c("tidyverse","BiocManager")); BiocManager::install(c("DESeq2","GEOquery","apeglm"))
# open rnaseq_de.R -> Session -> Set Working Directory -> To Source File Location -> Source

# Python (Colab or local Jupyter)
pip install pydeseq2 pandas matplotlib
# run rnaseq_de.ipynb
```

Data is fetched live from GEO (GSE147507). (GEO downloads are occasionally slow — rerun if a fetch times out.)

## 11. Abstract

*Performed differential-expression analysis of the SARS-CoV-2 airway host-response dataset (GSE147507) in
both R (DESeq2) and Python (PyDESeq2), fetched from the GEO API. Using a `~ cell_line + infection` design to
control for cell-line confounding, recovered ~3,976 differentially expressed genes dominated by inflammatory
cytokines and interferon-stimulated genes, and showed that controlling for the cell-line covariate increased
detected DE genes 4.5-fold — a concrete illustration of why the design formula determines correctness.
Cross-validated identical results across two independent DE toolchains.*

## Tech

`R` (DESeq2, GEOquery, apeglm) · `Python` (PyDESeq2, pandas) · GEO API · negative-binomial GLMs

## Files

```
rnaseq_de.R                     # R implementation (DESeq2)
rnaseq_de.ipynb                 # Python implementation (PyDESeq2)
results_R/deseq2_results.csv    # full DE results
results_R/volcano.png           # hero: shrunken-LFC volcano
results_R/pca.png               # sample PCA (separates by cell line)
results_py/…                    # Python equivalents
```
