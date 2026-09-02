# NeuroD2/NeuroD6 DKO — E15 cortex snRNA-seq

Analysis code for the supplemental single-nucleus RNA-seq in the Neurod2/Neurod6
manuscript.

Wild-type and Neurod2/Neurod6 double-knockout (DKO) E15 mouse cortex, one animal
per genotype. After QC filtering (>400 genes, <5% mitochondrial reads), doublet
removal with DoubletFinder, and removal of one low-quality cluster, the analysed
set is **12,882 WT and 26,385 DKO nuclei (39,267 total)**.

This repository contains **code only**. The input object and every generated
figure and table are distributed separately — see [Data](#data) and
[Outputs](#outputs).

---

## Contents

```
00_preprocessing_make_rds.R                            provenance for the input object
1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd  primary analysis  (res 0.5)
2_misspecified_clusters/scNeuroD_DKO_2misspecified.Rmd companion analysis (res 0.9)
```

### `00_preprocessing_make_rds.R`

How `scNeuroD_DKO_2025_doublet_Removed_azimuth.rds` was made: Cell Ranger
matrices → merge with genotype labels → QC metrics → filter → standard Seurat
workflow → per-genotype DoubletFinder → re-run workflow on the cleaned object
(this is the stored res-0.5 clustering) → Azimuth against `mousecortexref`.

It is **not** part of either notebook and is not re-run by them. It records
verbatim the code that produced the object, and is kept separate because it needs
things the analysis does not: the Cell Ranger outputs on a network volume, and
the `DoubletFinder`, `Azimuth`, `celldex` and `glmGamPoi` packages. A guard at
the top refuses to overwrite an existing object, so an accidental `source()`
cannot destroy the analysis input.

### The two analyses

Both notebooks read the **same object** and differ in exactly one thing: the
clustering resolution. Each writes only into its own directory.

| | `1_misspecified_cluster/` | `2_misspecified_clusters/` |
|---|---|---|
| Clustering | res 0.5, as stored in the object | re-clustered at res 0.9 |
| Named clusters | 20 | 21 |
| DKO-specific clusters | `Misspecified` (3,424 DKO / 1 WT) | `Misspecified UL` (3,522 / 1), `Misspecified DL` (2,109 / 96) |
| Role | **primary** — the resolution reported in the manuscript | companion — resolves the aberrant population in two |

They are deliberately **two files rather than one file behind a parameter**.
Each clustering needs its own hand-assigned numeric → cell-type mapping, and a
single file switching between them would silently mislabel every figure if
either clustering shifted. Instead each notebook *asserts* its own mapping
against known per-cluster cell counts in the `load_object` chunk, so a changed
clustering fails loudly rather than producing plausible, wrong output. The res-0.9
clustering and naming were verified to reproduce the archived 2026 run cell for
cell.

**Why both were run.** At res 0.5, Layer VIb appears expanded three-fold in the
mutant — 3.23% of WT nuclei against 10.04% of DKO nuclei. At res 0.9 that
cluster separates into a genuine Layer VIb, which is in fact *reduced*
(2.53% → 1.60%, 0.63×), and a DKO-specific Layer-VIb-like population enriched
10.7× (0.75% → 7.99%). Same nuclei, two descriptions. The apparent expansion is
real but is not an expansion of Layer VIb itself. The primary analysis is the
more conservative description — it does not subdivide a population on the basis
of one animal per genotype — and the companion analysis is what makes the nature
of that expansion interpretable.

---

## Data

The input object is **524 MB** and is not in this repository.

**Download:** [`scNeuroD_DKO_2025_doublet_Removed_azimuth.rds`](https://drive.google.com/file/d/1selJ4Y-o35d2nNGO4i4TDaafdXfyn5he/view?usp=sharing) (Google Drive)

From the command line, use `gdown` — a plain `curl` of the link returns Google's
virus-scan interstitial rather than the file, because it is over the 100 MB
scan limit:

```sh
pip install gdown
gdown 1selJ4Y-o35d2nNGO4i4TDaafdXfyn5he -O scNeuroD_DKO_2025_doublet_Removed_azimuth.rds
```

Verify the download before running anything — both notebooks assert this
checksum at the end of the run to confirm the object was not modified:

```
md5  be1d0e281dacfdfb2b26d655d8860250
```

```sh
md5 -q scNeuroD_DKO_2025_doublet_Removed_azimuth.rds   # macOS
md5sum scNeuroD_DKO_2025_doublet_Removed_azimuth.rds   # Linux
```

The object is **read-only** to both notebooks. Nothing in either file writes to
that path; derived objects are saved under new names.

---

## Running

Both notebooks build every path from two constants at the top of the file rather
than from `setwd()`, because knitr resets the working directory after every
chunk:

```r
PROJ_DIR     <- path.expand("~/SynologyDrive/NeuroD_2024/Revision July 28")
ANALYSIS_DIR <- file.path(PROJ_DIR, "1_misspecified_cluster")
```

`PROJ_DIR` is the path on the machine the archived results were produced on, and
has been **left unmodified** so the code matches the run that produced them.
Point it at your clone, put the `.rds` there, and the notebooks are
self-contained:

```r
PROJ_DIR <- path.expand("~/path/to/Neurod2d6DKO")
```

**Render with `Rscript`, not RStudio.** Forked parallelism requires a
non-interactive session; under RStudio the fork plan silently falls back to
`multisession`, which serialises the whole Seurat object to every worker.

```sh
Rscript -e 'rmarkdown::render("1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd")'
Rscript -e 'rmarkdown::render("2_misspecified_clusters/scNeuroD_DKO_2misspecified.Rmd")'
```

Each notebook writes `plots/` and `tables/` into its own directory and the HTML
report alongside the `.Rmd`. All are gitignored.

### Requirements

R with Seurat v5. Attached: `Seurat`, `ggplot2`, `dplyr`, `tidyr`, `tibble`,
`purrr`, `Matrix`, `ComplexHeatmap`, `circlize`, `future`, `future.apply`.

Called with `::` rather than attached, so they also work inside `multisession`
workers — which do not inherit attached packages: `clusterProfiler`,
`enrichplot`, `data.table`, `BiocParallel`, `ragg`, `scales`.

Must be installed but are never referenced directly: `org.Mm.eg.db` (passed to
`gseGO`/`enrichGO` as the *string* `"org.Mm.eg.db"`, so no SQLite connection is
ever serialised to a worker), `fgsea` (the estimator under `gseGO`), `magick`
(see rasterisation below), and `rmarkdown` / `knitr` to render.

`00_preprocessing_make_rds.R` additionally needs `DoubletFinder`
(github: chris-mcginnis-ucsf/DoubletFinder), `Azimuth`
(github: satijalab/azimuth), `celldex` and `glmGamPoi`.

### Machine notes

Written for 15 cores / 48 GB (10 heavy workers, 12 light, 16 GB
`future.globals.maxSize`); adjust `N_WORKERS_HEAVY` / `N_WORKERS_LIGHT` in the
config chunk. Stages touching the full Seurat object **fork**, so the object is
shared copy-on-write instead of serialised; stages calling AnnotationDbi get
fresh sessions instead, because forked processes inherit an open SQLite
connection that is not safe to use concurrently. Nested parallelism is disabled
inside workers — `fgsea` and `data.table` both grab every core by default, and
combined with an outer `future_lapply` that oversubscribes to 100+ threads and
runs slower than serial.

**Rasterisation.** `ComplexHeatmap` forces `type = "cairo"` whenever
`capabilities("cairo")` is `TRUE`; on macOS without XQuartz that flag is `TRUE`
but the DLL cannot load, and every rasterised heatmap dies. The notebooks probe
for a working device and prefer `ragg::agg_png`, which has no X11 dependency.
Install `magick` as well — without it `ComplexHeatmap` point-samples ~31,000 cell
columns into the raster image and most of the structure is thrown away.

---

## Statistics

**There is one animal per genotype.** This governs every choice below.

- **No p-values on genotype contrasts.** Cell-type proportions and pseudobulk
  expression differences are single observations, reported as descriptive effect
  sizes — proportions, ratios, log2 fold-changes. A cell-level test would treat
  nuclei from one animal as independent replicates, which is pseudoreplication,
  and would produce impressively small but meaningless p-values.
- **Pseudobulk (CPM) within-cluster DKO-vs-WT log2FC** for ranking and
  visualisation, in place of per-cell `FindMarkers` p-values.
- **One-vs-rest for genotype-private clusters.** The misspecified clusters have
  essentially no WT nuclei, so a within-cluster DKO-vs-WT contrast does not
  exist for them. They are compared against all other cells in the dataset — a
  marker-style contrast — consistently across DE, heatmaps, GSEA and ORA. This
  is a *different quantity* from the shared DKO-vs-WT columns and is captioned as
  such throughout.
- **Where a formal test is reported** (GSEA, GO over-representation), it tests
  gene-set structure within a ranked list, not a difference between animals.
- **GSEA ranks on per-cell `FindMarkers` log2FC**, which is continuous and so
  effectively tie-free; the p-values from that call are discarded and only the
  ordering is used. The tileplot *display* metric is the pseudobulk log2FC, which
  intentionally differs from the ranking metric — ranking needs a tie-free score,
  display needs consistency with the heatmaps.
- **Composition is computed within genotype**, so the ~2× difference in nuclei
  recovered between the two libraries cancels by construction.

---

## Outputs

Not in this repository — 3.8 GB, ~14,900 PDFs and 105 CSVs across the two
analyses, plus a rendered HTML report each (27 MB and 30 MB). Regenerate them by
rendering the notebooks, or request the archived copies.

```
plots/{dimplots,proportions,heatmaps,GSEA,enrichGO,Violin,gene_umaps,
       layer_markers,custom_genes,requested_markers,interneurons}/
tables/{celltype_proportions,cluster_cell_counts,Genotype_changes,QC_*,
        layer_*,interneuron_*,GSEA/,enrichGO/,interneurons/}
```
