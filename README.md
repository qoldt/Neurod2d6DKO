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

Three things live outside this repository: the raw sequencing reads and the Cell
Ranger outputs on NCBI GEO, and the analysis-ready Seurat object on Google Drive.
Only the last of these is needed to run the notebooks.

### Analysis-ready object (Google Drive)

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

### Raw and processed sequencing data (NCBI GEO)

<!-- GEO ACCESSIONS - fill in before release -->

Raw reads and Cell Ranger outputs are deposited in NCBI GEO under
**[GSEXXXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX)**,
BioProject **PRJNAXXXXXX**. Two samples, one animal per genotype, 10x Genomics
3' single-nucleus RNA-seq of E15 mouse cortex.

| GEO sample | Genotype | Nuclei after QC | Directory name the pipeline expects |
|---|---|---|---|
| GSMXXXXXXX | Wild type | 12,882 | `A4794_SP278_wt` |
| GSMXXXXXXX | Neurod2/Neurod6 DKO | 26,385 | `A4794_SP278_NeuroD2_D6DKO` |

**Raw FASTQ.** The 10x read files (`I1`, `R1`, `R2`) for each sample are in SRA,
linked from the corresponding GSM record and downloadable in bulk from the
BioProject:

```sh
prefetch SRRXXXXXXX
fasterq-dump --split-files SRRXXXXXXX
```

**Processed Cell Ranger output.** Each GSM carries the *filtered* feature-barcode
matrix as supplementary files — `barcodes.tsv.gz`, `features.tsv.gz`,
`matrix.mtx.gz` — produced by Cell Ranger `X.Y.Z` against the
`refdata-gex-mm10-XXXX-X` reference. These are the files
`00_preprocessing_make_rds.R` reads; the notebooks themselves never touch them.

**Rebuilding the object from GEO.** Arrange the downloaded matrices in the layout
Cell Ranger produces, since the script addresses them by sample directory:

```
$NEUROD_CELLRANGER_DIR/
├── A4794_SP278_wt/outs/filtered_feature_bc_matrix/{barcodes,features,matrix}
└── A4794_SP278_NeuroD2_D6DKO/outs/filtered_feature_bc_matrix/{barcodes,features,matrix}
```

```sh
export NEUROD_CELLRANGER_DIR=/path/to/cellranger_outputs
Rscript 00_preprocessing_make_rds.R
```

> **This will not reproduce the archived object bit for bit.** DoubletFinder's
> simulated doublets, UMAP and the Louvain clustering are all stochastic and the
> script sets no seed, so a rebuild yields a slightly different clustering. The
> notebooks assert exact per-cluster cell counts and will therefore stop rather
> than relabel — which is the intended behaviour. To reproduce the published
> figures, use the Google Drive object above and check its md5.

---

## Running

Clone it, put the `.rds` in the root, render. There are no absolute paths and
nothing to edit:

```sh
git clone https://github.com/qoldt/Neurod2d6DKO.git
cd Neurod2d6DKO
gdown 1selJ4Y-o35d2nNGO4i4TDaafdXfyn5he -O scNeuroD_DKO_2025_doublet_Removed_azimuth.rds

Rscript -e 'rmarkdown::render("1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd")'
Rscript -e 'rmarkdown::render("2_misspecified_clusters/scNeuroD_DKO_2misspecified.Rmd")'
```

Each notebook derives its paths from its own file location — `ANALYSIS_DIR` is
the directory holding the notebook, `PROJ_DIR` the repository root one level up.
They are constants rather than `setwd()` calls because knitr resets the working
directory after every chunk. Outputs stay inside the notebook's own directory;
the `.rds` is read from the root and never written to. If the object is missing,
the run stops immediately with download instructions instead of failing
somewhere downstream.

Two environment variables override the defaults; neither is normally needed:

| Variable | Effect |
|---|---|
| `NEUROD_PROJ_DIR` | Where the input object is looked for (default: repository root) |
| `NEUROD_RDS` | Full path to the object, to use a copy kept elsewhere |

**Render with `Rscript`, not RStudio** — as above. Forked parallelism requires a
non-interactive session; under RStudio the fork plan silently falls back to
`multisession`, which serialises the whole Seurat object to every worker.

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

`rmarkdown` needs **pandoc** on `PATH`. If you have RStudio but no standalone
pandoc, point at RStudio's bundled copy before rendering:

```sh
export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64
```

`00_preprocessing_make_rds.R` additionally needs `DoubletFinder`
(github: chris-mcginnis-ucsf/DoubletFinder), `Azimuth`
(github: satijalab/azimuth), `celldex` and `glmGamPoi`, plus
`NEUROD_CELLRANGER_DIR` pointing at the Cell Ranger run directory. It is
provenance, not part of either notebook, and refuses to overwrite an existing
object.

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

**Verified** on R 4.6.1 / Seurat 5.5.1, macOS aarch64, 15 cores / 48 GB: both
notebooks load the object and pass their cluster-naming assertions — 20 named
clusters at res 0.5, 21 at res 0.9 reproducing the archived run, 39,267 nuclei
in both.

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
