# Bulk RNA-seq Differential Expression — SARS-CoV-2 Host Response

Which genes respond to **SARS-CoV-2 infection** in human airway epithelial cells? A differential-expression
(DE) analysis of Blanco-Melo et al. 2020 (**GEO GSE147507**), implemented in **R (DESeq2)** and **Python
(PyDESeq2)** with matching results, then completed with **GO/KEGG functional enrichment (ORA + GSEA)** to turn
the gene list into pathways. The project doubles as a concrete lesson in why the **design formula** — not just
the code — decides whether a DE result is right.

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
7. **Functional enrichment**: GO/KEGG **over-representation** (ORA) on the up-regulated genes, and **GSEA** on
   all genes ranked by `sign(log2FC) × −log10(p)` — turning the gene list into pathways.

## 6. Results

| metric | value |
|--------|-------|
| significant genes, `~ cell_line + infection` (`padj<0.05`) | **3,976** |
| significant genes, `~ infection` only | 873 |
| effect of controlling for cell line | **4.5× more DE genes** |
| top up-regulated genes | IL36G, IL1A, SOD2, IRAK2, BIRC3, **MX1** |

Full results in [`results_R/deseq2_results.csv`](results_R/deseq2_results.csv). R and Python give identical
top genes and fold-changes (e.g. IL36G log₂FC ≈ 2.37, IL1A ≈ 3.25 in both).

**Enrichment.** ORA's top GO term is **cytokine-mediated signaling** (adjusted p ≈ 1×10⁻²⁶); GSEA ranks
**cytokine-cytokine receptor interaction, TNF, NF-κB, NOD-like receptor** and **JAK-STAT signalling** at the
top (NES ≈ +2.4–2.6, FDR q ≈ 0), while **oxidative phosphorylation** is the strongest down-regulated pathway
(NES ≈ −2.5). Tables in `enrichGO_up.csv` / `enrichKEGG_up.csv` / `gseKEGG.csv` (R) and `enrichr_up.csv` /
`gsea_prerank.csv` (Python).

## 7. Figures

`results_R/volcano.png` — shrunken-LFC volcano (hero, above); `pca.png` — sample PCA. Python equivalents in
`results_py/`.

## 8. Limitations

- Three cell lines with few replicates each — modest power; results are cell-culture, not in-vivo, biology.
- Genomic-span/count-level analysis only; no isoform- or allele-level resolution.
- Enrichment inherits annotation bias (well-studied immune genes are over-represented in GO/KEGG); ORA also
  depends on the significance threshold + chosen background, and pathways overlap (tests are not independent).
- DE + enrichment detect association with infection status, not causal/mechanistic direction.

## 9. Biological interpretation

The **PCA is dominated by cell line** (PC1 70% + PC2 24% ≈ 94% of variance); the infection effect is a subtle
within-cell-line shift. This is precisely why the **design formula** matters: with `~ infection` alone the
cell-line variance swamps the signal (873 genes), whereas `~ cell_line + infection` removes it and reveals
the infection response (**3,976 genes**). The top up-regulated genes are **inflammatory cytokines** (IL36G,
IL1A) and **antiviral interferon-stimulated genes** (MX1) — the imbalanced cytokine/interferon host response
to SARS-CoV-2 reported by Blanco-Melo et al.

**Functional enrichment closes the loop:** the gene list is not a random bag — it collapses into one coherent
program. ORA (top GO term *cytokine-mediated signaling*, adj p ≈ 1e-26) and GSEA (top pathways
*cytokine-cytokine receptor interaction, TNF, NF-κB, NOD-like receptor, JAK-STAT*, NES ≈ +2.5, FDR q ≈ 0) both
say infection drives a coordinated **innate-immune / inflammatory** response, while **oxidative
phosphorylation** is coordinately **repressed** (NES ≈ −2.5) — a known viral suppression of host energy
metabolism. That is the pathway-level conclusion the raw DEG list alone could not give.

## Language comparison (R vs Python)

Implemented in **DESeq2** (R) and **PyDESeq2** (Python) on the same data and design; the two agree to the
decimal on top genes and fold-changes — a strong cross-toolchain validation. The main idiomatic difference:
**DESeq2 expects counts as genes × samples, PyDESeq2 as samples × genes** (a transpose) — the most common
gotcha when porting a DE analysis between the two.

## 10. Reproduce

```bash
# R (RStudio)
install.packages(c("tidyverse","BiocManager"))
BiocManager::install(c("DESeq2","GEOquery","apeglm","clusterProfiler","org.Hs.eg.db","enrichplot"))
# open rnaseq_de.R -> Session -> Set Working Directory -> To Source File Location -> Source

# Python (Colab or local Jupyter)
pip install pydeseq2 gseapy pandas matplotlib
# run rnaseq_de.ipynb
```

Data is fetched live from GEO (GSE147507). (GEO downloads are occasionally slow — rerun if a fetch times out.)

## 11. Abstract

*Performed differential-expression analysis of the SARS-CoV-2 airway host-response dataset (GSE147507) in
both R (DESeq2) and Python (PyDESeq2), fetched from the GEO API. Using a `~ cell_line + infection` design to
control for cell-line confounding, recovered ~3,976 differentially expressed genes dominated by inflammatory
cytokines and interferon-stimulated genes, and showed that controlling for the cell-line covariate increased
detected DE genes 4.5-fold — a concrete illustration of why the design formula determines correctness. GO/KEGG
over-representation and GSEA resolved the gene list into a coordinated cytokine / NF-κB / JAK-STAT innate-immune
program with concomitant repression of oxidative phosphorylation. Cross-validated identical results across two
independent toolchains.*

## Tech

`R` (DESeq2, GEOquery, apeglm, clusterProfiler) · `Python` (PyDESeq2, gseapy, pandas) · GEO API ·
negative-binomial GLMs · GO/KEGG enrichment (ORA + GSEA)

## Files

```
rnaseq_de.R                     # R: DESeq2 DE + clusterProfiler enrichment (GO/KEGG ORA + GSEA)
rnaseq_de.ipynb                 # Python: PyDESeq2 DE + gseapy enrichment
results_R/deseq2_results.csv    # full DE results
results_R/volcano.png           # hero: shrunken-LFC volcano
results_R/pca.png               # sample PCA (separates by cell line)
results_R/enrichGO_up_dotplot.png, gseKEGG.csv, enrichGO_up.csv   # enrichment (R)
results_py/…                    # Python equivalents (volcano, enrichr_up.csv, gsea_prerank.csv)
```
