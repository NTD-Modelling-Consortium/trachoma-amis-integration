# Post-AMIS analysis — ground-truth plot export

Diagnostic plots comparing AMIS fits to survey data (trajectories, infections, parameter panels). Used for fit QA before downstream projection/scenario work.

## Scripts

| File | Role |
|------|------|
| [`export_ground_truth_plots.sh`](export_ground_truth_plots.sh) | Bash wrapper: validates inputs, stages legacy layout, runs R script, packages exports |
| [`plot_trachoma_fits_mtp.R`](plot_trachoma_fits_mtp.R) | R implementation (ggplot2 / patchwork); originally laptop-only, now env-driven |

Legacy exploratory notebooks/HTML live under [`summaries/`](summaries/) and are not part of the automated export path.

## When to run

After **fitting** completes. For plots that include projection-input checks, run after **projections-prep** as well (needs `InputBet_<IU>.csv` under `projections-prep/artefacts/projections/trachoma/<folder-id>/`).

Plot export is **not** invoked by `run_pipeline.py`. Run `export_ground_truth_plots.sh` explicitly once the pipeline artefacts you need are on disk (see [Running inside Docker](#running-inside-docker) below).

## Inputs

### From fitting-prep (usually image-baked)

- `Maps/trachoma_maps.rds`
- `Maps/trachoma_map_years.rds`
- `Maps/table_iu_idx_trachoma.csv`
- `trachoma/data/coverage/endgame_inputs/InputMDA_MTP_*.csv`

### From fitting (host-mounted for standalone plot export)

- `AMIS_output/amis_output<ID>.Rdata`
- `trajectories/trajectories_<ID>.Rdata`
- `infections/infections<ID>.Rdata`

These are written under `fitting/artefacts/` when the fitting stage runs. If you ran the pipeline in Docker with artefact directories mounted on the host (as `run_container.sh` does), point the wrapper at that same host tree.

### From projections-prep (host-mounted)

- `projections/trachoma/<folder-id>/<country>/<IU>/InputBet_<IU>.csv`

IUs with ESS below `--ess-threshold` (default 200) will not have InputBet files; exclude them via `FAILED_IDS` or accept missing batches.

## Outputs

Written to `PLOT_OUTPUT_DIR` (must be set explicitly):

| File | Description |
|------|-------------|
| `trachoma-ground-truth-plots-<RUN_DATE>.tar.gz` | Archive of everything under `post_AMIS_analysis/plots/` |
| `trachoma-ground-truth-plots-<RUN_DATE>.pdf` | Merged PDF (when the R script produced PDFs) |
| `trachoma-ground-truth-plots-<RUN_DATE>_preview.png` | First PNG copied when no PDFs (common in subset runs) |
| `trachoma-ground-truth-plots-<RUN_DATE>_5pct.pdf` | Optional; only if R output includes a `*5pct*` PDF |

Intermediate PNGs/PDFs live in `<workdir>/post_AMIS_analysis/plots/` during the run.

## Environment variables

### Wrapper (`export_ground_truth_plots.sh`)

See `export_ground_truth_plots.sh --help`. Key variables: `FOLDER_ID`, `PLOT_OUTPUT_DIR`, `FITTING_ARTEFACTS`, `FITTING_PREP_ARTEFACTS`, `PROJECTIONS_PREP_ARTEFACTS`, `INPUTBET_DIR_NAME`, `RUN_DATE`, `FAILED_IDS`, `KEEP_WORK_DIR`.

Pass **`PLOT_BATCH_IDS`** through to the R script (via `docker run -e` or the shell environment) when plotting a subset of batches.

### R script (`plot_trachoma_fits_mtp.R`)

| Variable | Description |
|----------|-------------|
| `TRACHOMA_PLOT_WORKDIR` | Root of staged workspace (wrapper sets this) |
| `TRACHOMA_INPUTBET_DIR` | Relative path to InputBet CSV folder under workdir |
| `PLOT_BATCH_IDS` | If set, plot only these TaskIDs and **exit before** aggregate map/beta sections |
| `FAILED_IDS` | TaskIDs to skip when `PLOT_BATCH_IDS` is unset |

## Running inside Docker

From the repository root, after pipeline stages have populated `fitting/artefacts/` and `projections-prep/artefacts/` on the host (for example via `./run_container.sh`):

```bash
FOLDER_ID=source-data-20250225-espen
RUN_DATE=20260819

docker run --rm \
  -e FOLDER_ID=${FOLDER_ID} \
  -e PLOT_OUTPUT_DIR=/plot-output \
  -e RUN_DATE=${RUN_DATE} \
  -e FITTING_ARTEFACTS=/ntdmc/trachoma-amis-integration/fitting/artefacts \
  -e PROJECTIONS_PREP_ARTEFACTS=/ntdmc/trachoma-amis-integration/projections-prep/artefacts \
  -v $(pwd)/fitting/artefacts:/ntdmc/trachoma-amis-integration/fitting/artefacts \
  -v $(pwd)/projections-prep/artefacts:/ntdmc/trachoma-amis-integration/projections-prep/artefacts \
  -v $(pwd)/plot-output:/plot-output \
  --entrypoint bash \
  trachoma-amis-pipeline:latest \
  -lc 'bash /ntdmc/trachoma-amis-integration/post_AMIS_analysis/export_ground_truth_plots.sh'
```

Add `-e PLOT_BATCH_IDS=10` to limit to one batch.

Create `plot-output/` on the host before running if it does not exist. Further detail: [`../README-docker-fitting.md`](../README-docker-fitting.md#ground-truth-plot-export-post-fit).

## Dependencies

Pre-installed in the Docker image: `reshape2`, `patchwork`, `viridis`, `weights`, plus the existing AMIS/plotting stack. The wrapper can install missing CRAN packages at runtime if needed.
