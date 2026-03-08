#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Read_Processing_parallel.sh - Parallel UMI extraction for smMIP processing
# Processes multiple samples concurrently using GNU parallel
# Memory-intensive R operations with controlled parallelism
# =============================================================================

# ----------- defaults -----------
BAM_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/bwa_primary"
PANEL_FILE="/Volumes/Seq_SSD/smMIP/10_23_25/Myeloid_Panel_Targets.chr.txt"
CODE_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/code/smMIP-tools-main/R"
OUT_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/ReadProcessing"
R_SCRIPT="${CODE_DIR}/map_smMIPs_extract_UMIs.R"

THREADS_PER_SAMPLE=2
PARALLEL_JOBS=4
OVERLAP=0.95
MAPQ=0
FORCE_REDO="${FORCE_REDO:-0}"

# Memory settings for parallel R sessions
# 30GB per instance × 4 = 120GB (leaves headroom from 128GB total)
R_MEM_PER_SAMPLE="30GB"

usage() {
  cat <<EOF
Usage:
  $0 --bam_dir DIR --panel FILE --code_dir DIR --out OUT_DIR [options]

Options:
  --bam_dir DIR        Input directory containing BAM files
  --panel FILE         Panel targets file
  --code_dir DIR       smMIP-tools R code directory
  --out OUT_DIR        Output directory
  --r_script FILE      Path to R script (default: CODE_DIR/map_smMIPs_extract_UMIs.R)
  --parallel N         Number of samples to process in parallel (default: 4)
  --threads N          Threads per sample for R (default: 2)
  --overlap N          Overlap threshold (default: 0.65)
  --mapq N             Minimum mapping quality (default: 0)
  --r_mem SIZE         R memory limit per sample (default: 30GB)
  --force              Force reprocessing even if outputs exist
  -h, --help           Show this help message

Resource allocation (default 4 parallel):
  - 2 threads/sample × 4 samples = 8 cores (R is memory-bound)
  - 30 GB RAM/sample × 4 samples = 120 GB RAM

Requires: GNU parallel (brew install parallel)

Note: R sessions are memory-intensive. If you experience memory issues,
reduce --parallel to 3 or increase --r_mem if you have more RAM.
EOF
}

say(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }

# ----------- args -----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam_dir) BAM_DIR="$2"; shift 2;;
    --panel) PANEL_FILE="$2"; shift 2;;
    --code_dir) CODE_DIR="$2"; R_SCRIPT="${CODE_DIR}/map_smMIPs_extract_UMIs.R"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --r_script) R_SCRIPT="$2"; shift 2;;
    --parallel) PARALLEL_JOBS="$2"; shift 2;;
    --threads) THREADS_PER_SAMPLE="$2"; shift 2;;
    --overlap) OVERLAP="$2"; shift 2;;
    --mapq) MAPQ="$2"; shift 2;;
    --r_mem) R_MEM_PER_SAMPLE="$2"; shift 2;;
    --force) FORCE_REDO=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -z "${BAM_DIR}" || -z "${PANEL_FILE}" || -z "${OUT_DIR}" ]] && { usage; exit 1; }

# ----------- checks -----------
need Rscript
need parallel

mkdir -p "${OUT_DIR}"
LOG_DIR="${OUT_DIR}/_logs"
mkdir -p "${LOG_DIR}"
PARALLEL_LOG="${LOG_DIR}/parallel_readproc.log"

# ----------- helper functions -----------
sample_done() {
  local s="$1" d="$2"
  [[ -s "${d}/${s}/${s}_clean.bam" ]] || return 1
  [[ -s "${d}/${s}/${s}_filtered_read_counts.txt" ]] || return 1
  [[ -s "${d}/${s}/${s}_raw_coverage_per_smMIP.txt" ]] || return 1
  [[ -s "${d}/${s}/${s}_UMI_usage_per_smMIP.txt" ]] || return 1
  return 0
}
export -f sample_done

# ----------- processing function -----------
process_one_sample() {
  local bam="$1"
  local base sample log status

  base="$(basename "$bam")"
  # Strip known BAM suffixes to derive sample name.
  # Handles names with hyphens (e.g., "KG-001" won't be truncated).
  # Try removing ".filtered.bam" first, then ".coorSorted.bam", then ".bam".
  sample="${base%.filtered.bam}"
  if [[ "$sample" == "$base" ]]; then
    sample="${base%.coorSorted.bam}"
  fi
  if [[ "$sample" == "$base" ]]; then
    sample="${base%.bam}"
  fi

  if [[ -z "$sample" || "$sample" == "$base" ]]; then
    echo "WARNING: Could not derive sample from: $base" >&2
    return 0
  fi

  # Skip if outputs exist (unless FORCE_REDO=1)
  if sample_done "$sample" "${_OUT_DIR}" && [[ "${_FORCE_REDO}" != "1" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} already processed - skipping."
    return 0
  fi

  mkdir -p "${_OUT_DIR}/${sample}"
  log="${_LOG_DIR}/${sample}.smmip.log"

  {
    echo "============================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing: ${sample}"
    echo "  BAM:  ${bam}"
    echo "  OUT:  ${_OUT_DIR}/${sample}"
    echo "  R_MAX_VSIZE=${_R_MEM} ; THREADS=${_THREADS} ; MAPQ=${_MAPQ} ; OVERLAP=${_OVERLAP}"
  } | tee "$log"

  # Run R with controlled memory
  R_MAX_VSIZE="${_R_MEM}" Rscript --no-restore --no-save \
    "${_R_SCRIPT}" \
    -b "$bam" \
    -p "${_PANEL_FILE}" \
    -s "$sample" \
    -o "${_OUT_DIR}" \
    -c "${_CODE_DIR}" \
    -t "${_THREADS}" \
    -O "${_OVERLAP}" \
    -M "${_MAPQ}" 2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}

  if [[ $status -ne 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAILED: ${sample} (exit $status). See: $log" >&2
    return $status
  fi

  if ! sample_done "$sample" "${_OUT_DIR}"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Expected outputs missing for ${sample}. See: $log" >&2
    return 3
  fi

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed: ${sample}"
}
export -f process_one_sample

# ----------- main -----------
say "Scanning for BAM files in ${BAM_DIR}..."

# Build list of BAM files
mapfile -d '' BAM_LIST < <(find "$BAM_DIR" -type f -name "*.filtered.bam" -print0)
[[ ${#BAM_LIST[@]} -gt 0 ]] || { echo "No .filtered.bam files found in ${BAM_DIR}"; exit 1; }

say "Found ${#BAM_LIST[@]} BAM files to process"
say "Running ${PARALLEL_JOBS} samples in parallel with ${THREADS_PER_SAMPLE} threads each"
say "R memory per sample: ${R_MEM_PER_SAMPLE}"
say "Total estimated memory: $((PARALLEL_JOBS * ${R_MEM_PER_SAMPLE%GB})) GB"

# Export variables for parallel subshells
export _OUT_DIR="${OUT_DIR}"
export _LOG_DIR="${LOG_DIR}"
export _THREADS="${THREADS_PER_SAMPLE}"
export _PANEL_FILE="${PANEL_FILE}"
export _CODE_DIR="${CODE_DIR}"
export _R_SCRIPT="${R_SCRIPT}"
export _OVERLAP="${OVERLAP}"
export _MAPQ="${MAPQ}"
export _R_MEM="${R_MEM_PER_SAMPLE}"
export _FORCE_REDO="${FORCE_REDO}"

# Run with GNU parallel
printf '%s\0' "${BAM_LIST[@]}" | \
  parallel -0 -j "${PARALLEL_JOBS}" \
    --bar \
    --joblog "${PARALLEL_LOG}" \
    --resume-failed \
    --halt soon,fail=1 \
    process_one_sample {}

say "All samples processed."
say "Job log: ${PARALLEL_LOG}"
