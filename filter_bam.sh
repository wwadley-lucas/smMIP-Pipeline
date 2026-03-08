#!/usr/bin/env bash
set -euo pipefail

# ----------- defaults -----------
BAM_DIR="./bwa_out"
OUT_DIR="./filtered_bam"
THREADS=6
MIN_MAPQ=50

usage() {
  cat <<EOF
Usage:
  $0 --bam_dir DIR --out OUT_DIR [--threads N] [--mapq N]

Filters BAM files:
  - Minimum mapping quality (default: 50)

EOF
}

say(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }

# ----------- args -----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam_dir) BAM_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mapq) MIN_MAPQ="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -z "${BAM_DIR}" || -z "${OUT_DIR}" ]] && { usage; exit 1; }

# ----------- checks -----------
need samtools
mkdir -p "${OUT_DIR}"

# ----------- main -----------
shopt -s nullglob
BAM_LIST=( "${BAM_DIR}"/*.bam )
[[ ${#BAM_LIST[@]} -gt 0 ]] || { echo "No .bam files found in ${BAM_DIR}"; exit 1; }

say "Filtering ${#BAM_LIST[@]} BAM file(s) with MAPQ >= ${MIN_MAPQ}"

for BAM in "${BAM_LIST[@]}"; do
  SAMPLE="$(basename "${BAM}" .bam)"
  OUT_BAM="${OUT_DIR}/${SAMPLE}.filtered.bam"

  say "Filtering ${SAMPLE}..."

  samtools view -b -q "${MIN_MAPQ}" -@ "${THREADS}" "${BAM}" -o "${OUT_BAM}"
  samtools index -@ "${THREADS}" "${OUT_BAM}"

  # Stats
  BEFORE=$(samtools view -c "${BAM}")
  AFTER=$(samtools view -c "${OUT_BAM}")
  PCT=$(awk -v after="$AFTER" -v before="$BEFORE" 'BEGIN {printf "%.1f", (after/before)*100}')

  say "  ${SAMPLE}: ${BEFORE} -> ${AFTER} reads (${PCT}% retained)"
done

say "All BAM files filtered. Output in ${OUT_DIR}"
