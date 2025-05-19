Scripts used for trachoma fitting and near term projections
================

Things to note: 

- Past version of projections was missing any MDA from 2022, and fits are missing many surveillance surveys. The latest version of these scripts (in `trachoma-amis-integration`) should be the correct version but they haven't been run
- The last runs in Feb 2025 we done on the cloud (<https://github.com/NTD-Modelling-Consortium/trachoma-docker-temp>) 
- May want to consider running with sigma=0.025 where there are IUs with ESS<200 or batches that failed (we have done this for other diseases)

## Scripts

- On HPC clusters, scripts that take long to run must be run through Slurm. Tables below 
show shell scripts because of this.

- Before running, clone <https://github.com/NTD-Modelling-Consortium/trachoma-amis-integration> and <https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma>.


### Preparing histories and maps for the fitting 

| R script  (see `Maps/` directory)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| prepare_histories_and_maps.R                                | runPrep.sh                    |
| prepare_histories_projections.R                             | runPrep_projections.sh        |


- **prepare_histories_and_maps.R**: produces the maps and the histories from 1996-2021.

- **prepare_histories_projections.R**: produces the histories from 1996-2025.

- Note: the histories (both for fitting and projections) are saved in the folder `ntd-model-trachoma/trachoma/data/coverage/endgame_inputs`. 


### Setup virtual environment

- The scripts were run from this Github <https://github.com/NTD-Modelling-Consortium/trachoma-amis-integration>. The model Github repository <https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma> must be installed separately since the remote repo does not have the MDA files required to run the fitting.

- Create a python virtual environment, clone and install both repositories (`trachoma-amis-integration` and `ntd-model-trachoma`). Note Feb 2025 runs were done using commits `main @ 5df8eb2` and `combineLatestChanges @ 4fda06b` respectively.

- Follow instructions in `trachoma-amis-integration` to install required R packages (doing this I still had issues with versions so for some packages I had to install the compatible version first manually and then I could install the remaining from `renv.lock`)


### Running the fitting


| R script  (see `Maps/` directory)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| trachoma_fitting.R                                | runFitting.sh                    |



### Plots for the model fits

- **plot_trachoma_fits_mtp.R**:  saves plots to summarise the fitting results and individual trajectories for each IU. Note I transferred the fitting results (`AMIS_output/`, `trajectories/`, `infections/` folders and the `InputBet_*` files) to my local machine and ran this, not on the cluster.


### Running the near term projections


- Before running these, we have to manually set the corresponding `failed_ids` in `preprocess_for_projections.R` and `realocate_InputPars_MTP.R` so we know which IUs failed to fit, and thus don't have any outputs to process. 

- In `realocate_InputPars_MTP.R` and `RunProjectionsTo2026.py` you should also update the folder name defined in the variable `folder_id` to reflect the current date.

- May also need to create directories `ntd-model-trachoma/projections` and `post_AMIS_analysis` if they don't already exist.

- Note that there is sometimes a maximum number of tasks that can be submitted at a time on HPC clusters (on the Warwick cluster it is 2000, I didn't reach the limit on the Oxford cluster). 


| R/Python script  (see `trachoma-amis-integration/` repo)    | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| preprocess_for_projections.R                                | samplePars.sh                 |         
| realocate_InputPars_MTP.R                                   | run_realocation_InputPars.sh  |
| RunProjectionsTo2026.py                                     | runProj.sh                    |
| RunProjectionsTo2026_stop_importation.py                     | NA                   |


- **preprocess_for_projections.R**: creates the 200 parameter vectors (simulated from the fitted 
models) used in projections.

- **realocate_InputPars_MTP.R**:  reorganises the files with the 200 samples used in projections, 
so that they are organised in the expected file hierarchy in the cloud.

- **RunProjectionsTo2026.py**: runs the near term projections for each IU.

- **RunProjectionsTo2026_stop_importation.py**: runs the near term projections for each IU but stops importation after last survey year


