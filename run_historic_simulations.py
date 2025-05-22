import argparse
import subprocess
from pathlib import Path
from typing import List
import os

import pandas as pd

# Script Behavior -
# 1. Load the table_iu_idx_trachoma.csv into DataFrame (Path = Maps)
# 2. Collect all the 0-indexed row-numbers matching the user specified task ID (identifies the IU codes)
# 3. Run the fitting and simulation for the given task ID and corresponding IUs

WORKING_DIR = Path(os.getenv("TRACHOMA_AMIS_DIR", ""))
PATH_TO_MAPS = WORKING_DIR / "Maps"
PATH_TO_TABLE_IU_IDX_TRACHOMA = PATH_TO_MAPS / "table_iu_idx_trachoma.csv"
PATH_TO_FITTING_SCRIPT = WORKING_DIR / "run_fit_and_process.sh"
PATH_TO_HISTORIC_SIMULATION_SCRIPT = WORKING_DIR / "run_proj.sh"
PATH_TO_INPUT_PARS_MTP_FILES = WORKING_DIR / "post_AMIS_analysis/InputPars_MTP_trachoma"


def run_script(path: Path, arg: int):
    command = [str(path), str(arg)]
    print(f"Running command: {' '.join(command)}")
    subprocess.run(command, check=True)
    print(f"Completed: {path} for argument {arg}")


def do_fitting(task_id: int):
    try:
        run_script(WORKING_DIR / PATH_TO_FITTING_SCRIPT, task_id)
        print(f"Fitting for task[{task_id}] completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error processing task[{task_id}]: {e}")


def do_historic_simulation(task_id: int, iu_indices: List[int]):
    for idx in iu_indices:
        try:
            run_script(WORKING_DIR / PATH_TO_HISTORIC_SIMULATION_SCRIPT, idx)
        except subprocess.CalledProcessError as e:
            print(f"Error processing task {task_id} with indices {iu_indices}: {e}")
    print(f"Task {task_id} with indices {iu_indices} completed successfully.")


def check_input_pars_files(iu_ids: List[int], indices: List[int]) -> set[int]:
    missing_indices = set()
    for iu, idx in zip(iu_ids, indices):
        input_file = (
            WORKING_DIR / PATH_TO_INPUT_PARS_MTP_FILES / f"InputPars_MTP_{iu}.csv"
        )
        if not input_file.exists():
            missing_indices.add(idx)

    return missing_indices


def main():
    expected_directory = WORKING_DIR
    current_directory = Path.cwd()
    if current_directory != expected_directory:
        raise RuntimeError(
            f"Error: This script must be run from '{expected_directory}', but it is running from '{current_directory}'."
        )

    parser = argparse.ArgumentParser(
        description="Run historic simulations and fitting."
    )
    parser.add_argument(
        "-t", "--task-id", type=int, help="Specify a TaskID to process", required=True
    )
    args = parser.parse_args()

    taskid_iu_df = pd.read_csv(
        PATH_TO_TABLE_IU_IDX_TRACHOMA, usecols=["IU_ID", "TaskID"]
    )

    if args.task_id:
        task_id = args.task_id
        task_rows = taskid_iu_df[taskid_iu_df["TaskID"] == task_id]
        iu_ids = task_rows["IU_ID"].tolist()
        indices = task_rows.index.tolist()
        do_fitting(task_id)

        # Due to the probabilistic nature of the fitting, some IUs may not have corresponding
        # parameter files, so we check here and exclude them.
        missing_indices = check_input_pars_files(iu_ids, indices)
        if missing_indices:
            print(
                f"Skipping IUs at {missing_indices} for task {task_id} (IU_IDS: {[iu_ids[indices.index(idx)] for idx in missing_indices]})."
            )
            indices = [idx for idx in indices if idx not in missing_indices]
        do_historic_simulation(task_id, indices)


if __name__ == "__main__":
    main()
