#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# filter_bam_parallel.sh - Parallel BAM filtering for smMIP processing
# Processes multiple samples concurrently using GNU parallel
# =============================================================================

# ----------- defaults -----------
BAM_DIR="./bwa_out"
OUT_DIR="./filtered_bam"
THREADS_PER_SAMPLE=6
PARALLEL_JOBS=4
MIN_MAPQ=50

usage() {
  cat <<EOF
Usage:
  $0 --bam_dir DIR --out OUT_DIR [options]

Options:
  --bam_dir DIR        Input directory containing BAM files
  --out OUT_DIR        Output directory for filtered BAM files
  --parallel N         Number of samples to process in parallel (default: 4)
  --threads N          Threads per sample for samtools (default: 6)
  --mapq N             Minimum mapping quality (default: 50)
  -h, --help           Show this help message

Filters BAM files:
  - Minimum mapping quality (default: 50)

Resource allocation (default 4 parallel):
  - 6 threads/sample × 4 samples = 24 cores
  - ~4 GB RAM/sample × 4 samples = ~16 GB RAM

Requires: GNU parallel (brew install parallel)
EOF
}

say(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }

# ----------- args -----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam_dir) BAM_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --parallel) PARALLEL_JOBS="$2"; shift 2;;
    --threads) THREADS_PER_SAMPLE="$2"; shift 2;;
    --mapq) MIN_MAPQ="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -z "${BAM_DIR}" || -z "${OUT_DIR}" ]] && { usage; exit 1; }

# ----------- checks -----------
need samtools
need parallel

mkdir -p "${OUT_DIR}"
LOG_DIR="${OUT_DIR}/_logs"
mkdir -p "${LOG_DIR}"
PARALLEL_LOG="${LOG_DIR}/parallel_filter.log"

# ----------- processing function -----------
filter_one_sample() {
  local BAM="$1"
  local SAMPLE OUT_BAM BEFORE AFTER PCT LOG

  SAMPLE="$(basename "${BAM}" .bam)"
  OUT_BAM="${_OUT_DIR}/${SAMPLE}.filtered.bam"
  LOG="${_LOG_DIR}/${SAMPLE}.filter.log"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtering ${SAMPLE}..." | tee "${LOG}"

  samtools view -b -q "${_MIN_MAPQ}" -@ "${_THREADS}" "${BAM}" -o "${OUT_BAM}" 2>> "${LOG}"
  samtools index -@ "${_THREADS}" "${OUT_BAM}" 2>> "${LOG}"

  # Stats
  BEFORE=$(samtools view -c "${BAM}")
  AFTER=$(samtools view -c "${OUT_BAM}")
  PCT=$(awk -v after="$AFTER" -v before="$BEFORE" 'BEGIN {printf "%.1f", (after/before)*100}')

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${SAMPLE}: ${BEFORE} -> ${AFTER} reads (${PCT}% retained)" | tee -a "${LOG}"
}
export -f filter_one_sample

# ----------- main -----------
shopt -s nullglob
BAM_LIST=( "${BAM_DIR}"/*.bam )
[[ ${#BAM_LIST[@]} -gt 0 ]] || { echo "No .bam files found in ${BAM_DIR}"; exit 1; }

say "Found ${#BAM_LIST[@]} BAM files to filter"
say "Running ${PARALLEL_JOBS} samples in parallel with ${THREADS_PER_SAMPLE} threads each"
say "Filtering with MAPQ >= ${MIN_MAPQ}"

# Export variables for parallel subshells
export _OUT_DIR="${OUT_DIR}"
export _LOG_DIR="${LOG_DIR}"
export _THREADS="${THREADS_PER_SAMPLE}"
export _MIN_MAPQ="${MIN_MAPQ}"

# Run with GNU parallel
printf '%s\0' "${BAM_LIST[@]}" | \
  parallel -0 -j "${PARALLEL_JOBS}" \
    --bar \
    --joblog "${PARALLEL_LOG}" \
    --resume-failed \
    --halt soon,fail=1 \
    filter_one_sample {}

say "All BAM files filtered. Output in ${OUT_DIR}"
say "Job log: ${PARALLEL_LOG}"
