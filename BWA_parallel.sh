#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# BWA_parallel.sh - Parallel BWA alignment for smMIP processing
# Processes multiple samples concurrently using GNU parallel
# =============================================================================

# ----------- defaults -----------
FASTQ_DIR="/Volumes/Seq_SSD/smMIP/KG001_01.22.25/RAW FASTQ copy"
REF_FA="/Volumes/Seq_SSD/smMIP/Universal/hg19/hg19.fa"
OUT_DIR="./bwa_out"
THREADS_PER_SAMPLE=6
PARALLEL_JOBS=4
PL="ILLUMINA"

R1_GLOB="*R1*.fastq.gz"
R2_TOKEN="R2"

BWA_BIN="${BWA_BIN:-bwa}"

usage() {
  cat <<EOF
Usage:
  $0 --fastq_dir DIR --ref REF.fa --out OUT_DIR [options]

Options:
  --fastq_dir DIR      Input directory containing FASTQ files
  --ref REF.fa         Reference genome FASTA file
  --out OUT_DIR        Output directory for BAM files
  --parallel N         Number of samples to process in parallel (default: 4)
  --threads N          Threads per sample for BWA/samtools (default: 6)
  --r1_glob GLOB       Glob pattern for R1 files (default: *R1*.fastq.gz)
  --r2_token TOKEN     Token to replace R1 with for R2 (default: R2)
  -h, --help           Show this help message

Resource allocation (default 4 parallel):
  - 6 threads/sample × 4 samples = 24 cores
  - ~8 GB RAM/sample × 4 samples = ~32 GB RAM

Requires: GNU parallel (brew install parallel)
EOF
}

say(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }

# ----------- args -----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq_dir) FASTQ_DIR="$2"; shift 2;;
    --ref) REF_FA="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --parallel) PARALLEL_JOBS="$2"; shift 2;;
    --threads) THREADS_PER_SAMPLE="$2"; shift 2;;
    --r1_glob) R1_GLOB="$2"; shift 2;;
    --r2_token) R2_TOKEN="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -z "${FASTQ_DIR}" || -z "${REF_FA}" || -z "${OUT_DIR}" ]] && { usage; exit 1; }

# ----------- checks -----------
need "${BWA_BIN}"
need samtools
need parallel

mkdir -p "${OUT_DIR}/_logs" "${OUT_DIR}/_tmp"
LOG_DIR="${OUT_DIR}/_logs"
PARALLEL_LOG="${LOG_DIR}/parallel_bwa.log"

# Index reference if needed
[[ -f "${REF_FA}.bwt" || -f "${REF_FA}.0123" ]] || { say "Indexing ${REF_FA} ..."; "${BWA_BIN}" index "${REF_FA}"; }
[[ -f "${REF_FA}.fai" ]] || samtools faidx "${REF_FA}"

# ----------- processing function -----------
process_one_sample() {
  local R1="$1"
  local base R2 SAMPLE RG COOR_BAM LOG SDIR

  base="$(basename "${R1}")"

  # Derive R2
  if [[ "${base}" == *R1* ]]; then
    R2="$(dirname "$R1")/${base/R1/${_R2_TOKEN}}"
  else
    R2="$(dirname "$R1")/${base/READ1/READ2}"
  fi

  if [[ ! -f "${R2}" ]]; then
    echo "WARNING: Missing R2 for ${base} (expected ${R2}); skipping."
    return 0
  fi

  # Sample name
  SAMPLE="${base}"
  SAMPLE="${SAMPLE/.fastq.gz/}"
  SAMPLE="${SAMPLE/.fq.gz/}"
  SAMPLE="${SAMPLE/.txt.gz/}"
  SAMPLE="${SAMPLE/.R1/}"
  SAMPLE="${SAMPLE/_R1/}"
  SAMPLE="${SAMPLE/READ1/}"

  SDIR="${_OUT_DIR}"
  LOG="${_LOG_DIR}/${SAMPLE}.bwa.log"
  : > "${LOG}"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning ${SAMPLE}" | tee -a "${LOG}"
  echo "  R1: ${R1}" | tee -a "${LOG}"
  echo "  R2: ${R2}" | tee -a "${LOG}"

  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${_PL}"
  COOR_BAM="${SDIR}/${SAMPLE}.bam"

  "${_BWA_BIN}" mem -M -Y -t "${_THREADS}" -R "${RG}" "${_REF_FA}" "${R1}" "${R2}" \
    2>> "${LOG}" \
  | samtools sort -@ "${_THREADS}" -o "${COOR_BAM}" -T "${_OUT_DIR}/_tmp/${SAMPLE}.tmp" -

  samtools index -@ "${_THREADS}" "${COOR_BAM}"
  samtools flagstat -@ "${_THREADS}" "${COOR_BAM}" > "${SDIR}/${SAMPLE}.flagstat.txt"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed ${SAMPLE} -> ${COOR_BAM}" | tee -a "${LOG}"
}
export -f process_one_sample

# ----------- main -----------
shopt -s nullglob
R1_LIST=( "${FASTQ_DIR}"/${R1_GLOB} )
[[ ${#R1_LIST[@]} -gt 0 ]] || { echo "No R1 files at ${FASTQ_DIR}/${R1_GLOB}"; exit 1; }

say "Found ${#R1_LIST[@]} R1 files to process"
say "Running ${PARALLEL_JOBS} samples in parallel with ${THREADS_PER_SAMPLE} threads each"
say "Total threads: $((PARALLEL_JOBS * THREADS_PER_SAMPLE))"

# Export variables for parallel subshells
export _BWA_BIN="${BWA_BIN}"
export _REF_FA="${REF_FA}"
export _OUT_DIR="${OUT_DIR}"
export _LOG_DIR="${LOG_DIR}"
export _THREADS="${THREADS_PER_SAMPLE}"
export _PL="${PL}"
export _R2_TOKEN="${R2_TOKEN}"

# Run with GNU parallel
printf '%s\0' "${R1_LIST[@]}" | \
  parallel -0 -j "${PARALLEL_JOBS}" \
    --bar \
    --joblog "${PARALLEL_LOG}" \
    --resume-failed \
    --halt soon,fail=1 \
    process_one_sample {}

say "All samples processed. Job log: ${PARALLEL_LOG}"
