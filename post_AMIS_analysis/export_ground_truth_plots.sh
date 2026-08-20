#!/usr/bin/env bash
# Stage fit pipeline artefacts into the legacy layout expected by plot_trachoma_fits_mtp.R,
# run the plotting script, and export canonical plot bundles.
set -euo pipefail

FOLDER_ID="${FOLDER_ID:-}"
PLOT_OUTPUT_DIR="${PLOT_OUTPUT_DIR:-}"
RUN_DATE="${RUN_DATE:-$(date +%Y%m%d)}"
FITTING_PREP_ARTEFACTS="${FITTING_PREP_ARTEFACTS:-}"
FITTING_ARTEFACTS="${FITTING_ARTEFACTS:-}"
PROJECTIONS_PREP_ARTEFACTS="${PROJECTIONS_PREP_ARTEFACTS:-}"
INPUTBET_DIR_NAME="${INPUTBET_DIR_NAME:-}"
WORK_DIR="${WORK_DIR:-}"
R_SCRIPT="${R_SCRIPT:-}"
FAILED_IDS="${FAILED_IDS:-}"
KEEP_WORK_DIR="${KEEP_WORK_DIR:-false}"

CONTAINER_AMIS_ROOT="${CONTAINER_AMIS_ROOT:-/ntdmc/trachoma-amis-integration}"
DEFAULT_R_SCRIPT="${CONTAINER_AMIS_ROOT}/post_AMIS_analysis/plot_trachoma_fits_mtp.R"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required (flag or environment variable):
  --folder-id ID                 FOLDER_ID / pipeline --folder-id
  --output-dir PATH              Directory for exported plot bundles
  --fitting-artefacts PATH       Host/container path to fitting/artefacts

Optional:
  --fitting-prep-artefacts PATH  Default: ${CONTAINER_AMIS_ROOT}/fitting-prep/artefacts
  --projections-prep-artefacts PATH
                                 Default: ${CONTAINER_AMIS_ROOT}/projections-prep/artefacts
  --inputbet-dir-name NAME       Subdir under post_AMIS_analysis/ for InputBet CSVs
                                 Default: <folder-id>-input-bet-csv
  --run-date YYYYMMDD            Used in export filenames (default: today)
  --work-dir PATH                Staging directory (default: temp dir, removed on exit)
  --r-script PATH                Plot script (default: ${DEFAULT_R_SCRIPT})
  --failed-ids LIST              Comma/space-separated TaskIDs to skip (passed to R)
  --keep-work-dir                Do not remove the staging directory on exit
  -h, --help

Environment mirrors the same names (FOLDER_ID, PLOT_OUTPUT_DIR, etc.).
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_dir() {
    local path="$1"
    local label="$2"
    [[ -d "${path}" ]] || die "missing ${label}: ${path}"
}

require_nonempty_dir() {
    local path="$1"
    local label="$2"
    require_dir "${path}" "${label}"
    if [[ -z "$(ls -A "${path}" 2>/dev/null || true)" ]]; then
        die "empty ${label}: ${path}"
    fi
}

require_glob() {
    local pattern="$1"
    local label="$2"
    # shellcheck disable=SC2086
    compgen -G "${pattern}" >/dev/null || die "missing ${label}: no match for ${pattern}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --folder-id)
            FOLDER_ID="$2"
            shift 2
            ;;
        --output-dir)
            PLOT_OUTPUT_DIR="$2"
            shift 2
            ;;
        --run-date)
            RUN_DATE="$2"
            shift 2
            ;;
        --fitting-prep-artefacts)
            FITTING_PREP_ARTEFACTS="$2"
            shift 2
            ;;
        --fitting-artefacts)
            FITTING_ARTEFACTS="$2"
            shift 2
            ;;
        --projections-prep-artefacts)
            PROJECTIONS_PREP_ARTEFACTS="$2"
            shift 2
            ;;
        --inputbet-dir-name)
            INPUTBET_DIR_NAME="$2"
            shift 2
            ;;
        --work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --r-script)
            R_SCRIPT="$2"
            shift 2
            ;;
        --failed-ids)
            FAILED_IDS="$2"
            shift 2
            ;;
        --keep-work-dir)
            KEEP_WORK_DIR=true
            shift
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1 (try --help)"
            ;;
    esac
done

[[ -n "${FOLDER_ID}" ]] || die "--folder-id / FOLDER_ID is required"
[[ -n "${PLOT_OUTPUT_DIR}" ]] || die "--output-dir / PLOT_OUTPUT_DIR is required"
[[ -n "${FITTING_ARTEFACTS}" ]] || die "--fitting-artefacts / FITTING_ARTEFACTS is required"

FITTING_PREP_ARTEFACTS="${FITTING_PREP_ARTEFACTS:-${CONTAINER_AMIS_ROOT}/fitting-prep/artefacts}"
PROJECTIONS_PREP_ARTEFACTS="${PROJECTIONS_PREP_ARTEFACTS:-${CONTAINER_AMIS_ROOT}/projections-prep/artefacts}"
INPUTBET_DIR_NAME="${INPUTBET_DIR_NAME:-${FOLDER_ID}-input-bet-csv}"
R_SCRIPT="${R_SCRIPT:-${DEFAULT_R_SCRIPT}}"

require_dir "${FITTING_PREP_ARTEFACTS}/Maps" "fitting-prep Maps"
require_file="${FITTING_PREP_ARTEFACTS}/Maps/trachoma_maps.rds"
[[ -f "${require_file}" ]] || die "missing fitting-prep map: ${require_file}"
[[ -f "${FITTING_PREP_ARTEFACTS}/Maps/trachoma_map_years.rds" ]] \
    || die "missing fitting-prep map years: ${FITTING_PREP_ARTEFACTS}/Maps/trachoma_map_years.rds"
[[ -f "${FITTING_PREP_ARTEFACTS}/Maps/table_iu_idx_trachoma.csv" ]] \
    || die "missing IU lookup: ${FITTING_PREP_ARTEFACTS}/Maps/table_iu_idx_trachoma.csv"

require_nonempty_dir "${FITTING_ARTEFACTS}/AMIS_output" "fitting AMIS_output"
require_nonempty_dir "${FITTING_ARTEFACTS}/trajectories" "fitting trajectories"
require_nonempty_dir "${FITTING_ARTEFACTS}/infections" "fitting infections"
require_glob "${FITTING_PREP_ARTEFACTS}/trachoma/data/coverage/endgame_inputs/InputMDA_MTP_*.csv" \
    "fitting-prep InputMDA_MTP CSVs"

inputbet_root="${PROJECTIONS_PREP_ARTEFACTS}/projections/trachoma/${FOLDER_ID}"
require_dir "${inputbet_root}" "projections-prep InputBet tree (${FOLDER_ID})"
inputbet_count="$(find "${inputbet_root}" -type f -name 'InputBet_*.csv' | wc -l | tr -d ' ')"
[[ "${inputbet_count}" -gt 0 ]] \
    || die "no InputBet_*.csv under ${inputbet_root}"

[[ -f "${R_SCRIPT}" ]] || die "plot script not found: ${R_SCRIPT}"

cleanup_work_dir() {
    if [[ "${KEEP_WORK_DIR}" == true ]]; then
        return 0
    fi
    if [[ -n "${WORK_DIR}" && "${WORK_DIR}" == /tmp/* && -d "${WORK_DIR}" ]]; then
        rm -rf "${WORK_DIR}"
    fi
}

if [[ -z "${WORK_DIR}" ]]; then
    WORK_DIR="$(mktemp -d)"
    trap cleanup_work_dir EXIT
else
    mkdir -p "${WORK_DIR}"
fi

echo "Staging plot workspace: ${WORK_DIR}"

stage_tree() {
    local src="$1"
    local dest="$2"
    mkdir -p "${dest}"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "${src}/" "${dest}/"
    else
        cp -a "${src}/." "${dest}/"
    fi
}

mkdir -p "${WORK_DIR}/Maps"
mkdir -p "${WORK_DIR}/AMIS_output"
mkdir -p "${WORK_DIR}/trajectories"
mkdir -p "${WORK_DIR}/infections"
mkdir -p "${WORK_DIR}/ntd-model-trachoma/trachoma/data/coverage/endgame_inputs"
mkdir -p "${WORK_DIR}/post_AMIS_analysis/plots"
inputbet_stage="${WORK_DIR}/post_AMIS_analysis/${INPUTBET_DIR_NAME}"
mkdir -p "${inputbet_stage}"

stage_tree "${FITTING_PREP_ARTEFACTS}/Maps" "${WORK_DIR}/Maps"
stage_tree "${FITTING_ARTEFACTS}/AMIS_output" "${WORK_DIR}/AMIS_output"
stage_tree "${FITTING_ARTEFACTS}/trajectories" "${WORK_DIR}/trajectories"
stage_tree "${FITTING_ARTEFACTS}/infections" "${WORK_DIR}/infections"
stage_tree \
    "${FITTING_PREP_ARTEFACTS}/trachoma/data/coverage/endgame_inputs" \
    "${WORK_DIR}/ntd-model-trachoma/trachoma/data/coverage/endgame_inputs"

while IFS= read -r -d '' bet_file; do
    cp -a "${bet_file}" "${inputbet_stage}/"
done < <(find "${inputbet_root}" -type f -name 'InputBet_*.csv' -print0)

echo "Staged ${inputbet_count} InputBet CSV(s) into post_AMIS_analysis/${INPUTBET_DIR_NAME}/"

export TRACHOMA_PLOT_WORKDIR="${WORK_DIR}"
export TRACHOMA_INPUTBET_DIR="post_AMIS_analysis/${INPUTBET_DIR_NAME}"
if [[ -n "${FAILED_IDS}" ]]; then
    export FAILED_IDS
fi

echo "Running plot script: ${R_SCRIPT}"
Rscript -e 'pkgs <- c("reshape2", "patchwork", "viridis"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) { cat("Error: missing R packages:", paste(missing, collapse=", "), "\nRebuild the trachoma-amis-pipeline Docker image.\n", file=stderr); quit(status=1) }'
Rscript "${R_SCRIPT}"

mkdir -p "${PLOT_OUTPUT_DIR}"
plots_dir="${WORK_DIR}/post_AMIS_analysis/plots"
tarball="${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}.tar.gz"
combined_pdf="${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}.pdf"

echo "Packaging plot outputs into ${PLOT_OUTPUT_DIR}"
tar -czf "${tarball}" -C "${plots_dir}" .

mapfile -t pdf_files < <(find "${plots_dir}" -maxdepth 1 -type f -name '*.pdf' | sort)
if [[ ${#pdf_files[@]} -eq 0 ]]; then
    mapfile -t pdf_files < <(find "${plots_dir}" -maxdepth 1 -type f -name '*.png' | sort)
fi
if [[ ${#pdf_files[@]} -eq 0 ]]; then
    die "plot script produced no PDF/PNG files under ${plots_dir}"
fi

if command -v pdfunite >/dev/null 2>&1 && [[ "${pdf_files[0]}" == *.pdf ]]; then
    pdfunite "${pdf_files[@]}" "${combined_pdf}"
elif command -v gs >/dev/null 2>&1 && [[ "${pdf_files[0]}" == *.pdf ]]; then
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${combined_pdf}" "${pdf_files[@]}"
elif [[ "${pdf_files[0]}" == *.png ]]; then
    cp -a "${pdf_files[0]}" "${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_preview.png"
    # No PDF merge for PNG-only smoke outputs; tarball still holds all plot files.
    combined_pdf=""
else
    cp -a "${pdf_files[0]}" "${combined_pdf}"
    echo "Warning: pdfunite/gs not found; copied ${pdf_files[0]} as ${combined_pdf}" >&2
fi

while IFS= read -r -d '' pct_pdf; do
    cp -a "${pct_pdf}" "${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_5pct.pdf"
    echo "Included 5pct PDF: $(basename "${pct_pdf}")"
    break
done < <(find "${plots_dir}" -maxdepth 1 -type f -iname '*5pct*.pdf' -print0)

echo "Exported:"
echo "  ${tarball}"
if [[ -n "${combined_pdf}" && -f "${combined_pdf}" ]]; then
    echo "  ${combined_pdf}"
fi
if [[ -f "${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_preview.png" ]]; then
    echo "  ${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_preview.png"
fi
if [[ -f "${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_5pct.pdf" ]]; then
    echo "  ${PLOT_OUTPUT_DIR}/trachoma-ground-truth-plots-${RUN_DATE}_5pct.pdf"
fi
