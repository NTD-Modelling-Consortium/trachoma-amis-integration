Trachoma AMIS Integration
=========================

Provides the relevant wrapper code to be able to call the trachoma
simulation from the AMIS loop to compute fitted parameters.

# Usage

## Prerequistes

 * R
 * Renv
 * CMake
 * Python >= 3.9

## Installation

Start with creating a Python virtual environment and activating it:

```shell
python3 -m venv trachoma-venv
source trachoma-venv/bin/activate
```

With the virtual environment activated, install the present AMIS
integration package:

```shell
git clone https://github.com/NTD-Modelling-Consortium/trachoma-amis-integration.git
cd trachoma-amis-integration/
pip install .
```

Finally, install R dependencies:


1. Ensure the dependencies are locked for the version of R you are
  using, by running the `setup_r_env.R` script: `Rscript setup_r_env.R`.
2. Launch a R session `R`
3. Create a R virtual environment with `renv::init()`
4. Install dependencies with `renv::restore()`

## Running

From the root of the project:

```bash
$ Rscript trachoma_fitting.R
```

This should sensibly pick the correct virtual environment if activated.
You can always specify the path to the virtual enviroment with:

```
export RETICULATE_PYTHON_ENV=/path/to/virtual/env/
```

## Upgrading AMIS

AMIS is locked to specific commit in the renv.lock file.

If the version on the repo is unchanged, then using `renv::update` won't work.

Instead, you must remove and re-add it:

```
> renv::remove("AMISforInfectiousDiseases")
> renv::install("evandrokonzen/AMISforInfectiousDiseases-dev")
> renv::snapshot()
```