#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Level_Bas_Calls_parallel.sh - Parallel pileup generation for smMIP processing
# Processes multiple samples concurrently using GNU parallel
# Memory-intensive R operations with controlled parallelism
# =============================================================================

# ----------- defaults -----------
READPROC_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/ReadProcessing"
PANEL_FILE="/Volumes/Seq_SSD/smMIP/10_23_25/Myeloid_Panel_Targets.chr.txt"
CODE_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/code/smMIP-tools-main/R"
R_SCRIPT="${CODE_DIR}/smMIP_level_raw_and_consensus_pileups.R"

OUT_PILEUP="${READPROC_DIR}/pileup"

THREADS_PER_SAMPLE=4
PARALLEL_JOBS=4

# Pileup settings
MMQ=20
MBQ=25
RANK="F"
UMI="T"
FAM_SIZE=2
CONS_CUTOFF=0.7

FORCE_REDO="${FORCE_REDO:-0}"

# Memory settings for parallel R sessions
# 25GB per instance × 4 = 100GB (leaves headroom from 128GB total)
R_MEM_PER_SAMPLE="25GB"

usage() {
  cat <<EOF
Usage:
  $0 --readproc_dir DIR --panel FILE --code_dir DIR [options]

Options:
  --readproc_dir DIR   Input directory from Read_Processing step
  --panel FILE         Panel targets file
  --code_dir DIR       smMIP-tools R code directory
  --out OUT_DIR        Output directory for pileups (default: READPROC_DIR/pileup)
  --r_script FILE      Path to R script (default: CODE_DIR/smMIP_level_raw_and_consensus_pileups.R)
  --parallel N         Number of samples to process in parallel (default: 4)
  --threads N          Threads per sample for R (default: 4)
  --mmq N              Minimum mapping quality (default: 20)
  --mbq N              Minimum base quality (default: 25)
  --rank F|T           Report top alt only if T (default: F)
  --umi T|F            Include UMI columns (default: T)
  --fam_size N         Family size threshold (default: 2)
  --cons_cutoff N      Consensus cutoff (default: 0.7)
  --r_mem SIZE         R memory limit per sample (default: 25GB)
  --force              Force reprocessing even if outputs exist
  -h, --help           Show this help message

Resource allocation (default 4 parallel):
  - 4 threads/sample × 4 samples = 16 cores
  - 25 GB RAM/sample × 4 samples = 100 GB RAM

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
    --readproc_dir) READPROC_DIR="$2"; OUT_PILEUP="${READPROC_DIR}/pileup"; shift 2;;
    --panel) PANEL_FILE="$2"; shift 2;;
    --code_dir) CODE_DIR="$2"; R_SCRIPT="${CODE_DIR}/smMIP_level_raw_and_consensus_pileups.R"; shift 2;;
    --out) OUT_PILEUP="$2"; shift 2;;
    --r_script) R_SCRIPT="$2"; shift 2;;
    --parallel) PARALLEL_JOBS="$2"; shift 2;;
    --threads) THREADS_PER_SAMPLE="$2"; shift 2;;
    --mmq) MMQ="$2"; shift 2;;
    --mbq) MBQ="$2"; shift 2;;
    --rank) RANK="$2"; shift 2;;
    --umi) UMI="$2"; shift 2;;
    --fam_size) FAM_SIZE="$2"; shift 2;;
    --cons_cutoff) CONS_CUTOFF="$2"; shift 2;;
    --r_mem) R_MEM_PER_SAMPLE="$2"; shift 2;;
    --force) FORCE_REDO=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -z "${READPROC_DIR}" || -z "${PANEL_FILE}" ]] && { usage; exit 1; }

# ----------- checks -----------
need Rscript
need parallel

mkdir -p "${OUT_PILEUP}"
LOG_DIR="${OUT_PILEUP}/_logs"
mkdir -p "${LOG_DIR}"
PARALLEL_LOG="${LOG_DIR}/parallel_pileup.log"

# ----------- helper functions -----------
sample_pileup_done() {
  local sample="$1"
  [[ -s "${_OUT_PILEUP}/${sample}_raw_pileup.txt" ]] && return 0
  [[ -s "${_OUT_PILEUP}/${sample}_consensus_pileup.txt" ]] && return 0
  return 1
}
export -f sample_pileup_done

# ----------- processing function -----------
process_one_sample() {
  local sample_dir="$1"
  local sample bam log status

  sample="$(basename "$sample_dir")"
  [[ -z "$sample" ]] && return 0

  # Prefer <SAMPLE>_clean.bam; otherwise use newest *clean.bam
  bam="${sample_dir%/}/${sample}_clean.bam"
  if [[ ! -f "$bam" ]]; then
    bam="$(find "$sample_dir" -type f -name "*clean.bam" -print0 2>/dev/null \
           | xargs -0 ls -1t 2>/dev/null | head -n 1 || true)"
  fi

  if [[ -z "${bam:-}" || ! -f "$bam" ]]; then
    echo "WARNING: No clean BAM found for sample '${sample}' under: $sample_dir" >&2
    return 0
  fi

  # Skip if outputs exist (unless FORCE_REDO=1)
  if sample_pileup_done "$sample" && [[ "${_FORCE_REDO}" != "1" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} already processed - skipping."
    return 0
  fi

  log="${_LOG_DIR}/${sample}.pileup.log"

  {
    echo "============================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing: ${sample}"
    echo "  BAM: $bam"
    echo "  OUT: ${_OUT_PILEUP}"
    echo "  R_MAX_VSIZE=${_R_MEM} ; THREADS=${_THREADS} ; MMQ=${_MMQ} ; MBQ=${_MBQ}"
    echo "  RANK=${_RANK} ; UMI=${_UMI} ; FAM_SIZE=${_FAM_SIZE} ; CONS_CUTOFF=${_CONS_CUTOFF}"
  } | tee "$log"

  # Run R with controlled memory
  R_MAX_VSIZE="${_R_MEM}" Rscript --no-restore --no-save "${_R_SCRIPT}" \
    -b "$bam" \
    -p "${_PANEL_FILE}" \
    -s "$sample" \
    -o "${_OUT_PILEUP}" \
    -c "${_CODE_DIR}" \
    -m "${_MMQ}" \
    -q "${_MBQ}" \
    -t "${_THREADS}" \
    -r "${_RANK}" \
    -u "${_UMI}" \
    -f "${_FAM_SIZE}" \
    -v "${_CONS_CUTOFF}" 2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}

  if [[ $status -ne 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAILED: ${sample} (exit $status). See: $log" >&2
    return $status
  fi

  if ! sample_pileup_done "$sample"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Expected pileup outputs not found for ${sample}. See: $log" >&2
    return 3
  fi

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed: ${sample}"
}
export -f process_one_sample

# ----------- main -----------
say "Scanning for sample directories in ${READPROC_DIR}..."

# Build list of sample directories
mapfile -d '' SAMPLE_DIRS < <(find "$READPROC_DIR" -mindepth 1 -maxdepth 1 -type d -print0)
[[ ${#SAMPLE_DIRS[@]} -gt 0 ]] || { echo "No sample directories found in ${READPROC_DIR}"; exit 1; }

say "Found ${#SAMPLE_DIRS[@]} sample directories to process"
say "Running ${PARALLEL_JOBS} samples in parallel with ${THREADS_PER_SAMPLE} threads each"
say "R memory per sample: ${R_MEM_PER_SAMPLE}"
say "Total estimated memory: $((PARALLEL_JOBS * ${R_MEM_PER_SAMPLE%GB})) GB"

# Export variables for parallel subshells
export _OUT_PILEUP="${OUT_PILEUP}"
export _LOG_DIR="${LOG_DIR}"
export _THREADS="${THREADS_PER_SAMPLE}"
export _PANEL_FILE="${PANEL_FILE}"
export _CODE_DIR="${CODE_DIR}"
export _R_SCRIPT="${R_SCRIPT}"
export _MMQ="${MMQ}"
export _MBQ="${MBQ}"
export _RANK="${RANK}"
export _UMI="${UMI}"
export _FAM_SIZE="${FAM_SIZE}"
export _CONS_CUTOFF="${CONS_CUTOFF}"
export _R_MEM="${R_MEM_PER_SAMPLE}"
export _FORCE_REDO="${FORCE_REDO}"

# Run with GNU parallel
printf '%s\0' "${SAMPLE_DIRS[@]}" | \
  parallel -0 -j "${PARALLEL_JOBS}" \
    --bar \
    --joblog "${PARALLEL_LOG}" \
    --resume-failed \
    --halt soon,fail=1 \
    process_one_sample {}

say "All samples processed."
say "Job log: ${PARALLEL_LOG}"
