# Trachoma Fitting and Near-Term Projections (Historical Simulations)
## Instructions for running in the Cloud using Docker

#### Things to note:
- Past version of projections was missing any MDA from 2022, and fits are missing many surveillance surveys. The latest version of these scripts (in `trachoma-amis-integration`) should be the correct version but they haven't been run
- The last runs in Feb 2025 were done on the cloud (<https://github.com/NTD-Modelling-Consortium/trachoma-docker-temp>) 
- May want to consider running with sigma=0.025 where there are IUs with ESS<200 or batches that failed (we have done this for other diseases) using `run_historic_simulations.sh --id=<failed_id> --folder_id="source-data-<yyyymmdd>" --amis_sigma=0.025 [--stop_importation]`.
- Keep note of which IUs had `ESS < 200` so that they can be excluded from the future projections (scenarios to 2040 run by Igor)

### Setup
- Clone this repo - [trachoma-amis-integration](https://github.com/NTD-Modelling-Consortium/trachoma-amis-integration).
- At the root of the repo, build the Docker image:
  ```shell
  DOCKER_BUILDKIT=1 docker build --ssh default=$SSH_AUTH_SOCK . -t trachoma-amis-pipeline
  ```
  This assumes that you have added the SSH keys on your system to your Github account.

### Pipeline
The pipeline can be considered to have three broad stages, each corresponding to some scripts -
1. **Preparation**:
    1. `Maps/prepare_histories_and_maps.R`: produces maps and histories from 1996-2021.
    2. `Maps/prepare_histories_projections.R`: produces the histories from 1996-2025.
2. **Fitting and Preprocessing Projections**:
    1. `trachoma_fitting.R`
    2. `preprocess_projections.R` - creates the `200` parameter vectors (simulated from the fitted models) used in projections.
    3. `realocate_InputPars_MTP.R` - reorganizes the files with the `200` samples used in projections, so that they are organized in the expected file hierarchy in the cloud.
3. **Near-term Projections**:
    1. `RunProjectionsTo2026.py` - runs the near-term projections for each IU, optionally with stopping importation after last survey year.

**NOTE**:
- The histories (for both fitting and projections) are saved inside the Docker container at `model/ntd-model-trachoma/trachoma/data/coverage/endgame_inputs`.
- All scripts are meant to be executed inside the Docker container.

The scripts in every stage can work with a single task/batch ID using the `--id` parameter. Some scripts expect the `SLURM_ARRAY_TASK_ID` environment variable to be defined if the `--id` argument is not provided. For convenience, a single entry point script `run_historical_simulations.sh` is provided which will execute the full pipeline. In normal usage, this is the only script that is needed.

Assuming invoking the script from inside the Docker container's shell, it's usage is as follows - 
```shell
./run_historic_simulations.sh --id=<id> --folder_id="source-data-<yyyymmdd>" [--stop_importation] [--amis_sigma=0.025]
```
This will produce the projections and place them inside the `model/ntd-model-trachoma/trachoma/projections`.

For more fine-grained usage where individual stages are invoked separately, pass the `--help` argument to the scripts to find out what arguments the script expects. This may be, particularly, useful when one or a subset of the stages is to be run. For example, refits might only need the *Preparation* and *Fitting* stages.

#### Docker
The pipeline can also be invoked using built container -
```shell
docker run -v ./trachoma:/ntdmc/trachoma-amis-integration/model/ntd-model-trachoma/projections/trachoma trachoma-amis-pipeline:latest --id=<id> --folder_id="source-data-<yyyymmdd>" [--stop_importation] [--amis_sigma=<0.025>]
```
This will produce the projections and place them inside the `trachoma` directory at the root of the repo on the host, allowing for easy access.
