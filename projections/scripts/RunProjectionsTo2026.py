#!/usr/bin/env python3
"""
Refactored script for running trachoma projections to 2026.

This script processes one or more Implementation Units (IUs) and generates
projections using fitted parameters from the AMIS process.
"""

import argparse
import logging
import os
import pickle
import sys
import time
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from ntdmc_trachoma.trachoma_functions import *

# Add trachoma_amis package to path
trachoma_amis_dir = os.getenv("TRACHOMA_AMIS_DIR", "")
if trachoma_amis_dir:
    sys.path.insert(0, trachoma_amis_dir)

from trachoma_amis.simulation_setup import (
    get_start_date,
    run_trachoma_simulations,
    PATH_TO_WORKING_DIR,
    PATH_TO_FITTING_PREP_ARTEFACTS,
    PATH_TO_PROJECTIONS_PREP_ARTEFACTS,
    PATH_TO_PROJECTIONS_ARTEFACTS,
)
from trachoma_amis.trachoma_params import (
    params as parameters,
    projection_params as sim_params,
    projection_config,
    demog,
)


def setup_logging(task_id: str) -> logging.Logger:
    """Set up logging for the projection run."""
    logging.basicConfig(
        level=logging.INFO,
        format=f"Task-{task_id}: %(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run trachoma projections to 2026",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-d", "--folder-id", required=True, help="Folder ID for outputs"
    )
    parser.add_argument(
        "-i",
        "--id",
        type=int,
        help="Row index corresponding to the SLURM array task ID",
    )
    parser.add_argument(
        "--num-cores", type=int, default=10, help="Number of CPU cores to use"
    )
    parser.add_argument(
        "--stop-importation",
        action="store_true",
        help="Stop importation of infections based on IU-specific year",
    )

    return parser.parse_args()


def validate_environment() -> None:
    """Validate that required environment variables and paths exist."""
    if not PATH_TO_WORKING_DIR.exists():
        raise ValueError(
            "TRACHOMA_AMIS_DIR environment variable not set or path doesn't exist"
        )


def load_iu_data(task_id: int) -> pd.DataFrame:
    """Load IU data for the given task ID."""
    iu_file = PATH_TO_FITTING_PREP_ARTEFACTS / "Maps" / "table_iu_idx_trachoma.csv"
    if not iu_file.exists():
        raise FileNotFoundError(f"IU index file not found: {iu_file}")

    df_iu_country = pd.read_csv(iu_file)
    batch_rows = df_iu_country[df_iu_country["TaskID"] == task_id]

    if batch_rows.empty:
        raise ValueError(f"No IUs found for task ID {task_id}")

    return batch_rows


def load_iu_parameters(iu_id: int, country: str, folder_id: str) -> dict:
    """Load fitted parameters for a specific IU."""
    param_file = (
        PATH_TO_PROJECTIONS_PREP_ARTEFACTS
        / "projections"
        / projection_config["species"]
        / folder_id
        / country
        / f"{country}{str(iu_id).zfill(5)}"
        / f"InputBet_{country}{str(iu_id).zfill(5)}.csv"
    )

    if not param_file.exists():
        raise FileNotFoundError(f"Parameter file not found: {param_file}")

    amisparams = pd.read_csv(param_file)
    amisparams.columns = [s.replace(" ", "") for s in amisparams.columns]

    # Extract parameter arrays
    seeds = amisparams.iloc[:, 0].astype(int).tolist()
    betas = amisparams.iloc[:, 1 : (int(sim_params["timesim"] / 52) + 1)].to_numpy()
    coverages = amisparams.iloc[:, -2].tolist()
    k_parameters = amisparams.iloc[:, -1].tolist()

    return {
        "seeds": seeds,
        "betas": betas,
        "coverages": coverages,
        "k_parameters": k_parameters,
    }


def calculate_stop_importation_time(stop_year: int, burnin: int) -> int:
    """Calculate simulation time to stop importation."""
    if stop_year is None:
        return -1

    stop_date = getOutputTimes([stop_year])
    return get_Intervention_times(stop_date, get_start_date(), burnin)[0]


def run_iu_projections(
    iu_id: int,
    country: str,
    folder_id: str,
    num_cores: int,
    stop_importation_year: Optional[int] = None,
) -> list:
    """Run projections for a single IU."""
    # Load parameters
    param_data = load_iu_parameters(iu_id, country, folder_id)

    # Set up MDA filepath
    mda_filepath = (
        PATH_TO_PROJECTIONS_PREP_ARTEFACTS
        / projection_config["species"]
        / "data"
        / "coverage"
        / "endgame_inputs"
        / f"InputMDA_MTP_projections_{iu_id}.csv"
    )

    # Calculate stop importation time
    timeToStopImportation = calculate_stop_importation_time(
        stop_importation_year, sim_params["burnin"]
    )

    # Convert betas to list format for the common function
    betas_list = [param_data["betas"][i, :] for i in range(len(param_data["seeds"]))]

    # Run simulations using common engine
    results = run_trachoma_simulations(
        seeds=param_data["seeds"],
        betas=betas_list,
        coverages=param_data["coverages"],
        k_parameters=param_data["k_parameters"],
        mda_filepath=mda_filepath,
        parameters=parameters,
        sim_params=sim_params,
        demog=demog,
        initial_infect_frac=projection_config["initial_infect_frac"],
        num_cores=num_cores,
        timeToStopImportation=timeToStopImportation,
        distToUse=projection_config["dist_to_use"],
    )

    return results


def getResultsNTDMC(results: list, start_date: date, burnin: int) -> pd.DataFrame:
    """Generate NTDMC results format from simulation results."""
    for i in range(len(results)):
        d = results[i][0].copy()
        prevs = np.array(d["True_Prev_Disease_children_1_9"])
        start = burnin
        step = 52
        chosenPrevs = prevs[start::step]

        if i == 0:
            df = pd.DataFrame(
                0, range(len(chosenPrevs)), columns=range(len(results) + 4)
            )
            df = df.rename(
                columns={0: "Time", 1: "age_start", 2: "age_end", 3: "measure"}
            )
            df.iloc[:, 0] = range(start_date.year, start_date.year + len(chosenPrevs))
            df.iloc[:, 1] = np.repeat(1, len(chosenPrevs))
            df.iloc[:, 2] = np.repeat(9, len(chosenPrevs))
            df.iloc[:, 3] = np.repeat("prevalence", len(chosenPrevs))

        df.iloc[:, i + 4] = chosenPrevs

    for i in range(len(results)):
        df = df.rename(columns={i + 4: "draw_" + str(i)})

    return df


def save_iu_results(results: list, iu_id: int, country: str, folder_id: str) -> None:
    """Save simulation results for an IU."""
    # Create output directory
    output_dir = (
        PATH_TO_PROJECTIONS_ARTEFACTS
        / projection_config["species"]
        / folder_id
        / country
        / f"{country}{str(iu_id).zfill(5)}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save simulation data as pickle
    sim_data = [item[0] for item in results]
    subset_sim_data = [
        {key: d[key] for key in projection_config["subset_keys"] if key in d}
        for d in sim_data
    ]

    pickle_file = (
        output_dir
        / f"{projection_config['species_prefix']}{country}{str(iu_id).zfill(5)}.p"
    )
    with open(pickle_file, "wb") as f:
        pickle.dump(subset_sim_data, f)

    # Save NTDMC format
    ntdmc_data = getResultsNTDMC(results, get_start_date(), sim_params["burnin"])
    csv_file = (
        output_dir
        / f"PrevDataset_{projection_config['species_prefix']}{country}{str(iu_id).zfill(5)}.csv"
    )
    ntdmc_data.to_csv(csv_file, index=False)


def process_single_iu(
    iu_data: pd.Series, folder_id: str, args: argparse.Namespace, logger: logging.Logger
) -> bool:
    """Process a single IU with error handling."""
    iu_id = iu_data["IU_ID"]
    country = iu_data["country"]

    try:
        logger.info(f"Starting projections for IU {iu_id} in {country}")

        # Get stop importation year if requested
        stop_year = None
        if args.stop_importation and "stop_importation_year" in iu_data:
            stop_year = iu_data["stop_importation_year"]

        # Run projections
        results = run_iu_projections(
            iu_id=iu_id,
            country=country,
            folder_id=folder_id,
            num_cores=args.num_cores,
            stop_importation_year=stop_year,
        )

        # Save results
        save_iu_results(results, iu_id, country, folder_id)

        logger.info(f"Completed projections for IU {iu_id}")
        return True

    except FileNotFoundError as e:
        logger.error(f"Missing input files for IU {iu_id}: {e}")
        return False
    except Exception as e:
        logger.error(f"Failed to process IU {iu_id}: {e}")
        return False


def main() -> int:
    """Main function for running projections."""
    start_time = time.time()

    try:
        # Parse arguments and validate environment
        args = parse_arguments()
        validate_environment()

        # Get task ID
        task_id = args.id if args.id else os.getenv("SLURM_ARRAY_TASK_ID")
        if not task_id:
            raise ValueError(
                "Must provide --id or set SLURM_ARRAY_TASK_ID environment variable"
            )
        task_id = int(task_id)

        # Set up logging
        logger = setup_logging(str(task_id))
        logger.info(f"Starting projections for task ID {task_id}")

        # Load IU data
        iu_data = load_iu_data(task_id)
        logger.info(f"Processing {len(iu_data)} IUs")

        # Process each IU
        success_count = 0
        for _, iu_row in iu_data.iterrows():
            if process_single_iu(iu_row, args.folder_id, args, logger):
                success_count += 1

        # Report results
        elapsed_time = time.time() - start_time
        logger.info(
            f"Completed: {success_count}/{len(iu_data)} IUs successful "
            f"in {elapsed_time:.1f} seconds"
        )

        return 0 if success_count == len(iu_data) else 1

    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
