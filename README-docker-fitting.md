# Trachoma Fitting and Near-Term Projections (Historical Simulations)
## Instructions for running in the Cloud using Docker

#### Things to note:
- Past version of projections was missing any MDA from 2022, and fits are missing many surveillance surveys. The latest version of these scripts (in `trachoma-amis-integration`) should be the correct version but they haven't been run
- The last runs in Feb 2025 were done on the cloud (<https://github.com/NTD-Modelling-Consortium/trachoma-docker-temp>) 
- May want to consider running with `amis-sigma=0.025` where there are IUs with `ESS < ess-threshold (default=200)` or batches that failed (we have done this for other diseases) using `bash run_container.sh --id=<failed_id> --folder-id="source-data-<yyyymmdd>" --amis-sigma=0.025 [--stop_importation]`.
- Keep note of which IUs had `ESS < ess-threshold (default=200)` so that they can be excluded from the future projections (scenarios to 2040 run by Igor)

### Setup
- Clone this repo - [trachoma-amis-integration](https://github.com/NTD-Modelling-Consortium/trachoma-amis-integration).
- At the root of the repo, build the Docker image:
  ```shell
  DOCKER_BUILDKIT=1 docker build --ssh default=$SSH_AUTH_SOCK . -t trachoma-amis-pipeline
  ```
  This assumes that you have added the SSH keys on your system to your Github account.

### Pipeline
The pipeline can be considered to have four broad stages, each corresponding to some scripts and generating specific output files:

#### 1. **Fitting Preparation** (`--stage=fitting-prep`)
**Script:** `prepare_histories_and_maps.R` - produces maps and histories from 1996-2021.

**Output Location:** `fitting-prep/artefacts/`

**Generated Files:**
- **Maps directory:**
  - `trachoma_espen_data.Rdata` - Processed ESPEN trachoma data
  - `trachoma_data.Rdata` - Aggregated trachoma data by IU/year
  - `iu_task_lookup.Rdata` - IU to TaskID lookup table
  - `trachoma_map_years.rds` - Years with survey data
  - `trachoma_maps.rds` - Processed trachoma maps with likelihood functions
  - `table_iu_idx_trachoma.csv` - IU/TaskID/country lookup table
- **Endgame inputs directory:**
  - `IUs_MTP_{id}.csv` - IU lists per batch
  - `InputMDA_MTP_{id}.csv` - MDA coverage data per batch

#### 2. **Fitting** (`--stage=fitting`)
**Script:** `trachoma_fitting.R` - runs AMIS fitting algorithm for each batch.

**Output Location:** `fitting/artefacts/`

**Generated Files:**
- `trajectories/trajectories_{id}.Rdata` - Simulated trajectories per batch
- `infections/infections{id}.Rdata` - Simulated infections per batch
- `AMIS_output/amis_output{id}.Rdata` - AMIS algorithm results per batch
- `ESS_NOT_REACHED.txt` - Log of batches where ESS target wasn't reached
- `summary.csv` - Summary statistics for all batch runs

#### 3. **Projections Preparation** (`--stage=projections-prep`)
**Scripts:**
- `prepare_histories_projections.R`: produces the histories from 1996-2025.
- `preprocess_for_projections.R`: creates the `amis-n-samples (default=200)` parameter vectors (simulated from the fitted models) used in projections.
- `realocate_InputPars_MTP.R`: reorganizes the files with the `amis-n-samples (default=200)` samples used in projections, so that they are organized in the expected file hierarchy in the cloud.

**Output Location:** `projections-prep/artefacts/`

**Generated Files:**
- `mda_history_trachoma.csv` - MDA history data
- `trachoma/data/coverage/endgame_inputs/IUs_MTP_{id}.csv` - IU lists per batch
- `trachoma/data/coverage/endgame_inputs/InputMDA_MTP_projections_{iu}.csv` - MDA coverage per IU
- `post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_{iu}.csv` - Parameters per IU
- `post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_{id}.rds` - Parameters per batch
- `post_AMIS_analysis/InputPars_MTP_trachoma/InputPars_MTP_allIUs.rds` - All IU parameters
- `projections/trachoma/{folder_id}/{country}/{country}{iu}/InputBet_{country}{iu}.csv` - Expanded parameters per IU

#### 4. **Near-term Projections** (`--stage=nearterm-projections`)
**Script:** `RunProjectionsTo2026.py` - runs the near-term projections for each IU, optionally with stopping importation after last survey year.

**Output Location:** `projections/artefacts/trachoma/{folder_id}/{country}/{country}{iu}/`

**Generated Files:**
- `Trachoma_{country}{iu}.p` - Simulation results (pickle format)
- `PrevDataset_Trachoma_{country}{iu}.csv` - Prevalence datasets (NTDMC format)

**NOTE**:
- The histories (for both fitting and projections) generated in the `fitting-prep` and `projections-prep` stages, are saved inside the Docker container at `fitting-prep/artefacts/trachoma/data/coverage/endgame_inputs` and `projections-prep/artefacts/trachoma/data/coverage/endgame_inputs` respectively.
- All scripts are meant to be executed inside the Docker container using the `run_container.sh` shell script.

#### Usage
The scripts in every stage can work with a single task/batch ID using the `--id` parameter. Some scripts expect the `SLURM_ARRAY_TASK_ID` environment variable to be defined if the `--id` argument is not provided. For convenience, a single entry point script `run_container.sh`, wrapping `run_pipeline.py`, is provided which will execute the full pipeline inside the Docker container and upon successful completion copy over the artefacts to the host. In normal usage, this is the only script that is needed.

It's usage is as follows - 
```shell
Usage: bash run_container.sh [options]

Required arguments:
  --id=<id>              SLURM batch/task ID for fitting and historic simulations
  --folder-id=<folder>   Folder for realocation (e.g., 'source-data-20250220')

Optional arguments:
  --stage                {fitting-prep, fitting, projections-prep, nearterm-projections, all, skip-fitting-prep} (default=all)
  --failed-ids=<ids>     Comma-separated list of failed batch/task IDs to skip
  --num-cores=<n>        Number of CPU cores to use for projections (default: 10)
  --stop-importation     Stop importation of infections based on IU-specific year
  --amis-sigma=<number>  AMIS 'sigma' parameter, expects a floating point number (default: 0.0025)
  --amis-target-ess      Target ESS parameter for AMIS (default: 500)
  --amis-n-samples       Number of AMIS samples (default: 1000)
  --ess-threshold        ESS threshold parameter (default: 200)
  --help                 Show this help message
```
For example,

```shell
bash run_container.sh --id=11 --folder-id="source-data-20250525" --num-cores=1 --stop-importation
```

This will produce the projections and place them inside the `projections/artefacts/trachoma`.

**NOTE**: `folder_id` is a bit of a misnomer because of the naming convention used for the directory `source-data-<yyyymmdd>`. The contents of this folder are outputs from the `fitting` and `projections-prep` stages but inputs to the subsequent `projections` stage.

For more fine-grained usage where individual stages are invoked separately, pass the `--help` argument to the scripts to find out what arguments are expected. This may be, particularly, useful when one or a subset of the stages is to be run. For example, refits might only need the `fitting-prep` and `fitting` stages.

#### Speeding up the pipeline (testing/development)
The AMIS related optional arguments can be specified to the `run_container.sh` script to speed up the pipeline during testing/development, as follows, for example - 

```shell
bash run_container.sh trachoma-amis-pipeline:latest \
  --id=<id> \
  --folder-id="source-data-<yyyymmdd>" \
  [--stop-importation] \
  [--amis-sigma=<0.025>] \
  [--amis-n-samples=10] \
  [--amis-target-ess=1] \
  [--ess-threshold=1]
```

