#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import enum
from pathlib import Path


class Stage(enum.Enum):
    FITTING_PREP = "fitting-prep"
    FITTING = "fitting"
    PROJECTIONS_PREP = "projections-prep"
    NEARTERM_PROJECTIONS = "nearterm-projections"
    ALL = "all"
    SKIP_FITTING_PREP = "skip-fitting-prep"


def validate_environment():
    """Validate required environment variables are set and we're in the correct directory."""
    trachoma_amis_dir = os.getenv("TRACHOMA_AMIS_DIR")
    if not trachoma_amis_dir:
        raise ValueError("TRACHOMA_AMIS_DIR environment variable is not set")

    current_dir = os.getcwd()
    if current_dir != trachoma_amis_dir:
        raise ValueError(
            f"This script must be run from '{trachoma_amis_dir}', "
            f"but it is running from '{current_dir}'."
        )


def run_command(command, description=None):
    """Run a command and handle any errors."""
    if description:
        print(f"{description}...")

    try:
        result = subprocess.run(command, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run the trachoma AMIS pipeline end-to-end"
    )

    # Required arguments
    parser.add_argument("--id", help="Batch/task ID to process")
    parser.add_argument(
        "--folder-id",
        required=True,
        type=str,
        help="Folder for realocation (e.g., 'source-data-20250220')",
    )

    # Optional arguments
    parser.add_argument(
        "--failed-ids",
        help="Comma-separated list of failed batch/task IDs to skip. Only used when --id is not specified.",
    )
    parser.add_argument(
        "--num-cores",
        type=int,
        default=10,
        help="Number of CPU cores to use for projections (default: 10)",
    )
    parser.add_argument(
        "--stop-importation",
        action="store_true",
        help="Stop importation of infections based on IU-specific year",
    )
    # AMIS-related parameters that can override environment variables
    parser.add_argument(
        "--amis-sigma",
        type=float,
        default=0.0025,
        help="AMIS 'sigma' parameter (default: 0.0025)",
    )
    parser.add_argument(
        "--amis-target-ess",
        type=int,
        help="Target ESS parameter for AMIS (default: 500)",
    )
    parser.add_argument(
        "--amis-n-samples",
        type=int,
        help="Number of AMIS samples (default: 1000)",
    )
    parser.add_argument(
        "--ess-threshold",
        type=int,
        help="ESS threshold parameter (default: 200)",
    )

    args = parser.parse_args()

    try:
        # Validate environment
        validate_environment()

        # Validate arguments
        if args.id and args.failed_ids:
            print("Warning: --failed-ids is ignored when --id is specified")

        r_args = ""
        # Export task ID for sub-processes if specified
        if args.id:
            os.environ["SLURM_ARRAY_TASK_ID"] = args.id
            r_args = f"--id={args.id}"

        # Step 1: Prepare histories and maps
        if not run_command(
            f"Rscript Maps/prepare_histories_and_maps.R {r_args}",
            "Preparing histories and maps",
        ):
            return 1

        if not run_command(
            f"Rscript Maps/prepare_histories_projections.R {r_args}",
            "Preparing histories projections",
        ):
            return 1

        # Step 2: Run fitting process
        # 2a: Trachoma fitting
        fit_args = []
        if args.id:
            fit_args.append(f"--id={args.id}")
        if args.amis_sigma:
            fit_args.append(f"--amis-sigma={args.amis_sigma}")
        if args.amis_target_ess:
            fit_args.append(f"--amis-target-ess={args.amis_target_ess}")
        if args.amis_n_samples:
            fit_args.append(f"--amis-n-samples={args.amis_n_samples}")
        if args.num_cores:
            fit_args.append(f"--num-cores={args.num_cores}")

        if not run_command(
            f"Rscript trachoma_fitting.R {' '.join(fit_args)}",
            "Running trachoma fitting",
        ):
            return 1

        # 2b: Preprocess projections
        preprocess_args = []
        if args.id:
            preprocess_args.append(f"--id={args.id}")
        elif args.failed_ids:
            preprocess_args.append(f"--failed-ids={args.failed_ids}")
        if args.ess_threshold:
            preprocess_args.append(f"--ess-threshold={args.ess_threshold}")

        if not run_command(
            f"Rscript preprocess_for_projections.R {' '.join(preprocess_args)}",
            "Preprocessing for projections",
        ):
            return 1

        # 2c: Realocate projections
        realocate_args = [f"--folder-id='{args.folder_id}'"]
        if args.id:
            realocate_args.append(f"--id={args.id}")
        elif args.failed_ids:
            realocate_args.append(f"--failed-ids={args.failed_ids}")
        if args.ess_threshold:
            realocate_args.append(f"--ess-threshold={args.ess_threshold}")

        if not run_command(
            f"Rscript realocate_InputPars_MTP.R {' '.join(realocate_args)}",
            "Realocating projections",
        ):
            return 1

        # /ntdmc/trachoma-amis-integration/model/ntd-model-trachoma/projections/trachoma/source-data-20250529/BFA/BFA05328/InputBet_BFA05328.csv
        PATH_TO_MODEL_DIR = Path(
            os.getenv("TRACHOMA_MODEL_DIR", "model/ntd-model-trachoma")
        )

        # Extract the countries corresponding to the `id` from the `Maps/table_iu_idx_trachoma.csv` file
        if args.id:
            import pandas as pd

            PATH_TO_WORKING_DIR = Path(os.getenv("TRACHOMA_AMIS_DIR", "."))
            PATH_TO_PROJECTIONS_PREP_ARTEFACTS = (
                PATH_TO_WORKING_DIR / "projections-prep" / "artefacts"
            )
            df_ius = pd.read_csv(
                PATH_TO_PROJECTIONS_PREP_ARTEFACTS / "table_iu_idx_trachoma.csv"
            )
            countries = df_ius[df_ius["TaskID"] == int(args.id)]["country"]
            iu_ids = df_ius[df_ius["TaskID"] == int(args.id)]["IU_ID"]
            for country, iu_id in zip(countries, iu_ids):
                input_bet_file = (
                    PATH_TO_PROJECTIONS_PREP_ARTEFACTS
                    / "projections"
                    / "trachoma"
                    / args.folder_id
                    / f"{country}/{iu_id}/InputBet_{iu_id}.csv"
                )
                if not input_bet_file.exists():
                    print(
                        f"InputBet file does not exist for {iu_id}. Skipping projections for this ID = {args.id}."
                    )
                    exit(1)

        # Step 3: Run historic simulations
        proj_args = [f"--folder-id='{args.folder_id}'"]
        if args.id:
            proj_args.append(f"--id={args.id}")
        if args.num_cores:
            proj_args.append(f"--num-cores={args.num_cores}")
        if args.stop_importation:
            proj_args.append("--stop-importation")

        if not run_command(
            f"python RunProjectionsTo2026.py {' '.join(proj_args)}",
            "Running historic simulations",
        ):
            return 1

        print("Pipeline completed successfully!")
        return 0

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
