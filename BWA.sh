#!/usr/bin/env bash
set -euo pipefail

# ----------- defaults -----------
FASTQ_DIR="/Users/lucaswadley/Desktop/RAW FASTQ copy"

# Reference (paper used hg19/GRCh37)
REF_FA="/Volumes/Seq_SSD/smMIP/Universal/hg19/hg19.fa"
OUT_DIR="./bwa_out"
THREADS=6
PL="ILLUMINA"

# File matching: adjust to your processed FASTQ names
R1_GLOB="*R1*.fastq.gz"   # e.g. SAMPLE.R1.trim.fq.gz, SAMPLE_R1.fastq.gz
R2_TOKEN="R2"             # how to derive R2 from R1 (replace the first 'R1' with 'R2')

BWA_BIN="${BWA_BIN:-bwa}" # use bwa-mem2 on x86 systems for 2-3x speedup

usage() {
  cat <<EOF
Usage:
  $0 --fastq_dir DIR --ref REF.fa --out OUT_DIR [--threads N]
     [--r1_glob GLOB] [--r2_token TOKEN]

Notes:
  - Assumes final, already-processed FASTQs (no trimming/UMI steps here).
  - Produces coordinate-sorted BAM (+ .bai).
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
    --threads) THREADS="$2"; shift 2;;
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
mkdir -p "${OUT_DIR}/_logs" "${OUT_DIR}/_tmp"
LOG_DIR="${OUT_DIR}/_logs"

# index reference if needed
[[ -f "${REF_FA}.bwt" || -f "${REF_FA}.0123" ]] || { say "Indexing ${REF_FA} ..."; "${BWA_BIN}" index "${REF_FA}"; }
[[ -f "${REF_FA}.fai" ]] || samtools faidx "${REF_FA}"

# ----------- main -----------
shopt -s nullglob
R1_LIST=( "${FASTQ_DIR}"/${R1_GLOB} )
[[ ${#R1_LIST[@]} -gt 0 ]] || { echo "No R1 files at ${FASTQ_DIR}/${R1_GLOB}"; exit 1; }

for R1 in "${R1_LIST[@]}"; do
  base="$(basename "${R1}")"
  # derive R2 by replacing the first occurrence of R1 with R2 token
  if [[ "${base}" == *R1* ]]; then
    R2="$(dirname "$R1")/${base/R1/${R2_TOKEN}}"
  else
    # fallback for READ1/READ2 style
    R2="$(dirname "$R1")/${base/READ1/READ2}"
  fi

  if [[ ! -f "${R2}" ]]; then
    echo "WARNING: Missing R2 for ${base} (expected ${R2}); skipping."
    continue
  fi

  # sample name from basename (strip common read token + extension variants)
  SAMPLE="${base}"
  SAMPLE="${SAMPLE/.fastq.gz/}"
  SAMPLE="${SAMPLE/.fq.gz/}"
  SAMPLE="${SAMPLE/.txt.gz/}"
  SAMPLE="${SAMPLE/.R1/}"
  SAMPLE="${SAMPLE/_R1/}"
  SAMPLE="${SAMPLE/READ1/}"

  SDIR="${OUT_DIR}"
  mkdir -p "${SDIR}"
  LOG="${LOG_DIR}/${SAMPLE}.bwa.log"; : > "${LOG}"

  say "▶ Aligning ${SAMPLE}" | tee -a "${LOG}"
  say "R1: ${R1}" | tee -a "${LOG}"
  say "R2: ${R2}" | tee -a "${LOG}"

  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${PL}"
  COOR_BAM="${SDIR}/${SAMPLE}.bam"

  # Skip if BAM and index already exist
  if [[ -f "${COOR_BAM}" && -f "${COOR_BAM}.bai" ]]; then
    say "Skipping ${SAMPLE} (BAM already exists: ${COOR_BAM})"
    echo "----------------------------------------------"
    continue
  fi

  "${BWA_BIN}" mem -M -Y -t "${THREADS}" -R "${RG}" "${REF_FA}" "${R1}" "${R2}" \
    2>> "${LOG}" \
  | samtools sort -@ "${THREADS}" -o "${COOR_BAM}" -T "${OUT_DIR}/_tmp/${SAMPLE}.tmp" -
  samtools index -@ "${THREADS}" "${COOR_BAM}"

  samtools flagstat -@ "${THREADS}" "${COOR_BAM}" > "${SDIR}/${SAMPLE}.flagstat.txt"
  say "✓ ${SAMPLE} → ${COOR_BAM}" | tee -a "${LOG}"
  echo "----------------------------------------------"
done

say "All samples processed."
