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


def get_args_for_fitting_prep(args):
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")

    return r_args


def get_args_for_fitting(args):
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")

    if args.amis_sigma:
        r_args.append(f"--amis-sigma={args.amis_sigma}")

    if args.amis_n_samples:
        r_args.append(f"--amis-n-samples={args.amis_n_samples}")

    if args.amis_target_ess:
        r_args.append(f"--amis-target-ess={args.amis_target_ess}")

    if args.num_cores:
        r_args.append(f"--num-cores={args.num_cores}")

    return r_args


def get_args_for_projections_prep(args):
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")
    elif args.failed_ids:
        r_args.append(f"--failed-ids={args.failed_ids}")

    if args.folder_id:
        r_args.append(f"--folder-id='{args.folder_id}'")

    if args.species:
        r_args.append(f"--species={args.species}")

    if args.ess_threshold:
        r_args.append(f"--ess-threshold={args.ess_threshold}")

    return r_args


def validate_fitting_prep_args(args):
    """Validate arguments required for fitting-prep stage."""
    missing = []
    if not args.id:
        missing.append("--id")

    if missing:
        raise ValueError(f"fitting-prep stage requires: {', '.join(missing)}")


def validate_fitting_args(args):
    """Validate arguments required for fitting stage."""
    missing = []
    if not args.id:
        missing.append("--id")

    if missing:
        raise ValueError(f"fitting stage requires: {', '.join(missing)}")


def validate_projections_prep_args(args):
    """Validate arguments required for projections-prep stage."""
    missing = []

    # Either --id or --failed-ids must be provided
    if not args.id and not args.failed_ids:
        missing.append("--id or --failed-ids")

    # --folder-id is required for realocate_InputPars_MTP.R script
    if not args.folder_id:
        missing.append("--folder-id")

    if missing:
        raise ValueError(f"projections-prep stage requires: {', '.join(missing)}")


def validate_nearterm_projections_args(args):
    """Validate arguments required for nearterm-projections stage."""
    missing = []
    if not args.id:
        missing.append("--id")
    if not args.folder_id:
        missing.append("--folder-id")

    if missing:
        raise ValueError(f"nearterm-projections stage requires: {', '.join(missing)}")


def get_args_for_nearterm_projections(args):
    r_args = []
    if args.folder_id:
        r_args.append(f"--folder-id='{args.folder_id}'")

    if args.id:
        r_args.append(f"--id={args.id}")

    if args.num_cores:
        r_args.append(f"--num-cores={args.num_cores}")

    if args.stop_importation:
        r_args.append("--stop-importation")

    return r_args


def do_fitting_prep(args):
    """Run the fitting preparation step of the pipeline."""
    validate_fitting_prep_args(args)

    PATH_TO_SCRIPTS = Path(
        os.getenv("PATH_TO_FITTING_PREP_SCRIPTS", "./fitting-prep/scripts")
    )

    r_args = get_args_for_fitting_prep(args)
    
    # Run prepare_histories_and_maps.R
    command = (
        f"Rscript {PATH_TO_SCRIPTS}/prepare_histories_and_maps.R {' '.join(r_args)}"
    )
    if not run_command(command, "Preparing histories and maps"):
        return False
    
    # Run prepare_histories_projections.R
    command = (
        f"Rscript {PATH_TO_SCRIPTS}/prepare_histories_projections.R {' '.join(r_args)}"
    )
    return run_command(command, "Preparing histories for projections")


def do_fitting(args):
    """Run the fitting step of the pipeline."""
    validate_fitting_args(args)

    PATH_TO_SCRIPTS = Path(os.getenv("PATH_TO_FITTING_SCRIPTS", "./fitting/scripts"))

    r_args = get_args_for_fitting(args)
    command = f"Rscript {PATH_TO_SCRIPTS}/trachoma_fitting.R {' '.join(r_args)}"
    return run_command(command, "Running trachoma fitting")


def do_projections_prep(args):
    """Run the projections preparation step of the pipeline."""
    validate_projections_prep_args(args)

    PATH_TO_SCRIPTS = Path(
        os.getenv("PATH_TO_PROJECTIONS_PREP_SCRIPTS", "./projections-prep/scripts")
    )

    scripts = [
        {
            "path": PATH_TO_SCRIPTS / "preprocess_for_projections.R",
            "args": set(["--id", "--ess-threshold", "--failed-ids", "--species"]),
        },
        {
            "path": PATH_TO_SCRIPTS / "realocate_InputPars_MTP.R",
            "args": set(
                [
                    "--id",
                    "--species",
                    "--folder-id",
                    "--ess-threshold",
                    "--failed-ids",
                ]
            ),
        },
    ]

    r_args = get_args_for_projections_prep(args)
    for script in scripts:
        script_args: set = script["args"]
        # Filter r_args to only include those that match script_args
        filtered_args = [
            arg
            for arg in r_args
            if any(arg.startswith(script_arg) for script_arg in script_args)
        ]
        print(f"Script arguments: {script_args}")
        print(f"Provided arguments: {r_args}")
        print(f"Filtered arguments for {script['path'].name}: {filtered_args}")
        command = f"Rscript {script['path']} {' '.join(filtered_args)}"
        if not run_command(command, f"Running {script['path'].name}"):
            return False
    return True


def do_nearterm_projections(args):
    """Run the near-term projections (to 2026) stage of the pipeline."""
    validate_nearterm_projections_args(args)

    PATH_TO_PROJECTIONS_SCRIPTS = Path(
        os.getenv("PATH_TO_PROJECTIONS_SCRIPTS", "./projections/scripts")
    )

    # Extract the countries corresponding to the `id` from the `Maps/table_iu_idx_trachoma.csv` file
    if args.id:
        import pandas as pd

        PATH_TO_WORKING_DIR = Path(os.getenv("TRACHOMA_AMIS_DIR", "."))
        PATH_TO_PROJECTIONS_PREP_ARTEFACTS = (
            PATH_TO_WORKING_DIR / "projections-prep" / "artefacts"
        )
        PATH_TO_FITTING_PREP_ARTEFACTS = (
            PATH_TO_WORKING_DIR / "fitting-prep" / "artefacts"
        )
        df_ius = pd.read_csv(
            PATH_TO_FITTING_PREP_ARTEFACTS / "Maps" / "table_iu_idx_trachoma.csv"
        )
        countries = df_ius[df_ius["TaskID"] == int(args.id)]["country"]
        iu_ids = df_ius[df_ius["TaskID"] == int(args.id)]["IU_ID"]
        for country, iu_id in zip(countries, iu_ids):
            iu_code = f"{country}{str(iu_id).zfill(5)}"
            input_bet_file = (
                PATH_TO_PROJECTIONS_PREP_ARTEFACTS
                / "projections"
                / "trachoma"
                / args.folder_id
                / f"{country}"
                / f"{iu_code}"
                / f"InputBet_{iu_code}.csv"
            )
            if not input_bet_file.exists():
                print(
                    f"InputBet file does not exist for {iu_code}. Skipping projections for this ID = {args.id}."
                )
                return False

    r_args = get_args_for_nearterm_projections(args)
    command = f"python {PATH_TO_PROJECTIONS_SCRIPTS}/RunProjectionsTo2026.py {' '.join(r_args)}"
    return run_command(command, "Running near-term projections to 2026")


STAGE_SEQUENCE_MAP = {
    Stage.FITTING_PREP.value: (do_fitting_prep,),
    Stage.FITTING.value: (do_fitting,),
    Stage.PROJECTIONS_PREP.value: (do_projections_prep,),
    Stage.NEARTERM_PROJECTIONS.value: (do_nearterm_projections,),
    Stage.ALL.value: (
        do_fitting_prep,
        do_fitting,
        do_projections_prep,
        do_nearterm_projections,
    ),
    Stage.SKIP_FITTING_PREP.value: (
        do_fitting,
        do_projections_prep,
        do_nearterm_projections,
    ),
}


def main():
    parser = argparse.ArgumentParser(
        description="Run the trachoma AMIS pipeline end-to-end",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
STAGE-SPECIFIC ARGUMENT REQUIREMENTS:

fitting-prep:
  Required: --id
  Optional: None
  Description: Prepares histories and maps for fitting stage

fitting:
  Required: --id
  Optional: --amis-sigma, --amis-n-samples, --amis-target-ess, --num-cores
  Description: Runs AMIS fitting process

projections-prep:
  Required: (--id OR --failed-ids) AND --folder-id
  Optional: --species, --ess-threshold
  Description: Prepares data for projections stage
  
nearterm-projections:
  Required: --id, --folder-id
  Optional: --num-cores, --stop-importation
  Description: Runs near-term projections to 2026

EXAMPLES:

  # Run individual stages
  %(prog)s --stage fitting-prep --id 1
  %(prog)s --stage fitting --id 1 --amis-sigma 0.003
  %(prog)s --stage projections-prep --id 1 --folder-id source-data-20250220
  %(prog)s --stage nearterm-projections --id 1 --folder-id source-data-20250220

  # Run full pipeline
  %(prog)s --id 1 --folder-id source-data-20250220

  # Skip fitting-prep (use existing fitting-prep artefacts)
  %(prog)s --stage skip-fitting-prep --id 1 --folder-id source-data-20250220

  # Process all IDs skipping the failed ones in projections-prep
  %(prog)s --stage projections-prep --failed-ids 1,2,3 --folder-id source-data-20250220
        """,
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--id",
        type=int,
        required=True,
        help="Batch/task ID to process",
    )

    parser.add_argument(
        "--folder-id",
        required=False,
        type=str,
        help="Folder for realocation (e.g., 'source-data-20250220'). Required for projections-prep and nearterm-projections stages.",
    )

    # Optional arguments
    parser.add_argument(
        "--failed-ids",
        type=str,
        required=False,
        help="Comma-separated list of failed batch/task IDs to skip. Only used when --id is not specified.",
    )

    parser.add_argument(
        "--num-cores",
        type=int,
        required=False,
        default=10,
        help="Number of CPU cores to use for projections (default: 10)",
    )

    parser.add_argument(
        "--stop-importation",
        required=False,
        action="store_true",
        help="Stop importation of infections based on IU-specific year",
    )
    # AMIS-related parameters that can override environment variables
    parser.add_argument(
        "--amis-sigma",
        type=float,
        required=False,
        default=0.0025,
        help="AMIS 'sigma' parameter (default: 0.0025)",
    )
    parser.add_argument(
        "--amis-target-ess",
        type=int,
        required=False,
        default=500,
        help="Target ESS parameter for AMIS (default: 500)",
    )
    parser.add_argument(
        "--amis-n-samples",
        type=int,
        required=False,
        default=1000,
        help="Number of AMIS samples (default: 1000)",
    )
    parser.add_argument(
        "--ess-threshold",
        type=int,
        required=False,
        default=200,
        help="ESS threshold parameter (default: 200)",
    )

    parser.add_argument(
        "--species",
        type=str,
        required=False,
        default="trachoma",
        help="Species to process (default: 'trachoma')",
    )
    parser.add_argument(
        "--stage",
        type=str,
        choices=[s.value for s in Stage],
        required=False,
        default=Stage.ALL.value,
        help="Stage of the pipeline to run. "
        f"Options: {', '.join(s.value for s in Stage)} (default: {Stage.ALL.value}). "
        "Each stage has different argument requirements - see STAGE-SPECIFIC ARGUMENT REQUIREMENTS below for details.",
    )

    args = parser.parse_args()

    try:
        # Validate environment
        validate_environment()

        # Set up environment if ID is provided
        if args.id:
            print(f"Running pipeline for ID: {args.id}")
            os.environ["SLURM_ARRAY_TASK_ID"] = str(args.id)
            print(f"SLURM_ARRAY_TASK_ID set to {os.environ['SLURM_ARRAY_TASK_ID']}")
            if args.failed_ids:
                print("Warning: --failed-ids is ignored when --id is specified")

        for stage_function in STAGE_SEQUENCE_MAP[args.stage]:
            print(f"Running stage: {stage_function.__name__}")
            try:
                if not stage_function(args):
                    print(f"Stage {stage_function.__name__} failed.", file=sys.stderr)
                    return 1
            except ValueError as e:
                print(
                    f"Stage {stage_function.__name__} argument validation failed: {e}",
                    file=sys.stderr,
                )
                return 1

        print("Pipeline completed successfully!")
        return 0
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
