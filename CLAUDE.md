# CLAUDE.md

Guidance for working in this repository. The README describes *what* the
analysis is; this file describes what will break if you are careless.

This is analysis code backing a manuscript revision. Its output goes into a
paper. A silently wrong figure is far worse than a loud failure — several of the
guards below exist because that already happened once.

## The one fact that governs everything

**There is one animal per genotype.** n = 1 WT, n = 1 DKO.

Never add a p-value to a between-genotype comparison. Not to cell-type
proportions, not to pseudobulk fold-changes, not "just for reference". A
cell-level test treats nuclei from a single animal as independent replicates —
pseudoreplication — and produces tiny, meaningless p-values that a reader will
take at face value. Genotype effects are reported as descriptive effect sizes
only: proportions, ratios, log2 fold-changes.

`FindMarkers` *is* called in the GSEA ranking chunks, but only its `avg_log2FC`
column is used, because the ranking needs a continuous tie-free score. The
p-values from those calls are discarded. Do not surface them.

GSEA and GO over-representation p-values are fine — they test gene-set structure
within a ranked list, not a difference between animals. Keep that distinction
intact in any caption you write.

## Two notebooks, deliberately duplicated

`1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd` (res 0.5, primary) and
`2_misspecified_clusters/scNeuroD_DKO_2misspecified.Rmd` (res 0.9, companion)
are ~2,300 lines each and differ in about 220 lines — the cluster name map,
`PRIVATE_CLUSTERS`, `CLUSTER_ORDER`, `EXPECTED_COUNTS`, and paths.

**Do not refactor them into one parameterised file.** This was considered and
rejected. Each clustering needs its own hand-assigned numeric → cell-type
mapping, and a single file switching between them would mislabel every figure
silently if either clustering shifted. The duplication is the safety mechanism.

**Any methodology change must be applied to both files.** They are supposed to
stay in lockstep on everything except the clustering. After editing one, diff
them and confirm the only differences are the intended ones plus the known
per-resolution config:

```sh
diff 1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd \
     2_misspecified_clusters/scNeuroD_DKO_2misspecified.Rmd
```

## Guards you must not "fix"

### The cluster-naming assertion

`load_object` maps cluster numbers to cell-type names by hand, then asserts the
resulting per-cluster cell counts against `EXPECTED_COUNTS` and stops the run on
any mismatch.

If that assertion fails, **the clustering has changed and every downstream label
is wrong.** Do not update `EXPECTED_COUNTS` to match the new numbers. That
converts a loud, correct failure into a plausible, wrong figure. Find out why
the clustering moved — wrong `.rds`, different Seurat version, changed
`DROP_CLUSTERS` — and fix that.

The same chunk also asserts that the misspecified cluster is genotype-private
(`WT <= 1`, `DKO > 3000`). The one-vs-rest contrast used throughout is only the
right comparison because that holds.

### The input object is read-only

`RDS_PATH` is only ever read. Nothing may `saveRDS()` to it. `session_info`
checks its md5 is still `be1d0e281dacfdfb2b26d655d8860250` at the end of the
run. Derived objects get new names.

`00_preprocessing_make_rds.R` has a guard refusing to overwrite an existing
object. Leave it. Do not run that script to "regenerate" anything — it needs
Cell Ranger outputs on a network volume and packages that are not installed
here.

### The per-cluster detection guard

The heatmap gene selection uses a cell-level detection-specificity guard,
`|pct_DKO − pct_WT| > 0.3`, computed **within each cluster**.

The predecessor computed it pooled across all cells. Pooling averages a
cluster-specific detection difference over the whole dataset and destroys it:
the largest pooled difference for any gene in this data is 0.2115, so the gate
admitted zero genes and every shared-cluster gene vanished from the heatmaps
without an error. Computed per cluster, the same threshold admits 132 genes.

The DE contrast is within-cluster. The guard must be too. Don't pool it.

## Shared vs private clusters

`SHARED_CLUSTERS` have both genotypes → within-cluster DKO-vs-WT pseudobulk.

`PRIVATE_CLUSTERS` (the misspecified ones) have essentially no WT nuclei, so a
within-cluster genotype contrast **does not exist**. They use one-vs-rest against
all other cells — a marker-style contrast — consistently across DE, heatmaps,
GSEA and ORA.

These are **different quantities in the same table and the same figure**. Any
caption, table header or summary you write must say so. Silently presenting a
one-vs-rest log2FC next to a DKO-vs-WT log2FC as if they were comparable is the
single easiest way to mislead a reader here.

Derive cluster scopes from `PRIVATE_CLUSTERS` / `CLUSTER_ORDER`, never by
hardcoding names. The predecessor hardcoded `"Misspecified"` in one place while
`case_when()` emitted `"Misspecified UL"` / `"Misspecified DL"`, and those bars
silently dropped out as `NA`.

## Conventions

- **Every figure goes through `save_plot()` / `save_pdf()`; every table through
  `save_table()`.** They guarantee the directory exists, the path is absolute,
  and nothing depends on ggplot2's `last_plot()` side effect. Never call
  `ggsave()` or `write.csv()` directly.
- **Paths come from `PROJ_DIR` / `ANALYSIS_DIR`, never `setwd()`.** knitr resets
  the working directory after every chunk.
- `PROJ_DIR` still points at the original SynologyDrive path. This is
  intentional — the code matches the run that produced the archived figures.
  Repoint it locally to run; don't commit a change to it without saying so.
- Chunks are quiet by default. A full run emits ~500 lines of incidental
  Seurat/ggplot2/clusterProfiler notes that bury the real output. Chunks whose
  log is worth reading opt in with `message=TRUE` in the header.
- Config lives in **one** place, the `config` chunk. The predecessor redefined
  genotype labels and cluster scopes in five chunks, which is how the
  `Misspecified` name drifted out of sync.

## Running

**Render with `Rscript`, never RStudio.** Forking requires a non-interactive
session; under RStudio `plan_fork()` silently falls back to `multisession` and
serialises the whole Seurat object to every worker.

```sh
Rscript -e 'rmarkdown::render("1_misspecified_cluster/scNeuroD_DKO_1misspecified.Rmd")'
```

A full run is long — thousands of PDFs, the `Violin` chunk alone dominates the
runtime. Don't launch one casually to check a small edit.

Parallel discipline, if you touch it:

- `plan_fork()` — stages touching the full Seurat object. Copy-on-write instead
  of serialising.
- `plan_session()` — stages calling AnnotationDbi/`org.Mm.eg.db`. Forked
  processes inherit an open SQLite connection that is not safe to use
  concurrently. `OrgDb` is passed as the *string* `"org.Mm.eg.db"` so no
  connection is ever serialised.
- `worker_setup()` — call it inside any new worker. `fgsea` and `data.table`
  each grab every core by default; combined with an outer `future_lapply` that
  oversubscribes to 100+ threads and runs slower than serial.
- Invert loops that would subset the Seurat object per iteration. The violin
  chunks fan out over *clusters* and loop genes inside — 10 subsets instead of
  ~5,000.

Heatmaps need a working raster device. `ComplexHeatmap` forces `type = "cairo"`
whenever `capabilities("cairo")` is `TRUE`; on macOS without XQuartz that flag
is `TRUE` but the DLL will not load and every rasterised heatmap dies. The
config chunk probes for a device and prefers `ragg::agg_png`. Keep `magick`
installed — without it ~31,000 cell columns are point-sampled into the raster
and most of the structure is thrown away.

## Repository hygiene

Code only. `.gitignore` covers `*.rds`, `plots/`, `tables/`, `*.html`. Do not
commit generated output — it is ~3.8 GB and ~14,900 PDFs — and do not commit the
524 MB input object; it is distributed via Google Drive.

## Before you claim something works

Rendering succeeds does not mean the analysis is right — the bugs that mattered
most here all produced clean runs with empty or mislabelled figures. If you
change gene selection, check how many genes each stage actually admits. If you
change cluster handling, check the counts. Report what you verified and what you
did not.
