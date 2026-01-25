#!/usr/bin/env bash
set -euo pipefail

# ── User paths ────────────────────────────────────────────────────────────────
BAM_DIR="/Volumes/Seq_SSD/smMIP/KG001_01.22.25/filtered_bam"
PANEL_FILE="/Volumes/Seq_SSD/smMIP/Universal/Myeloid_Panel_Targets.chr.txt"
CODE_DIR="/Volumes/Seq_SSD/smMIP/Universal/Code/R"
OUT_DIR="/Volumes/Seq_SSD/smMIP/KG001_01.22.25/Read_Processing"
R_SCRIPT="${CODE_DIR}/map_smMIPs_extract_UMIs.R"

THREADS=2
OVERLAP=0.95     # -O
MAPQ=50           # -M
FORCE_REDO="${FORCE_REDO:-0}"  # set to 1 to re-run even if outputs exist

mkdir -p "${OUT_DIR}"

LOG_DIR="${OUT_DIR}/_logs"
mkdir -p "${LOG_DIR}"

# ── Helpers ──────────────────────────────────────────────────────────────────
ensure_RMAX_in_Renviron() {
  local renv="$HOME/.Renviron"
  touch "$renv"
  if grep -q '^R_MAX_VSIZE=' "$renv"; then
    sed -i.bak '/^R_MAX_VSIZE=/d' "$renv"
  fi
  echo 'R_MAX_VSIZE=120GB' >> "$renv"
  export R_MAX_VSIZE="120GB"
  echo "ℹ️  R_MAX_VSIZE set to ${R_MAX_VSIZE}"
}

clear_R_environment() {
  R_MAX_VSIZE=120GB Rscript --no-restore --no-save -e "rm(list=ls()); invisible(gc());"
}

sample_done() {
  # $1 = sample id, $2 = outdir
  local s="$1" d="$2"
  [[ -s "${d}/${s}/${s}_clean.bam" ]]                         || return 1
  [[ -s "${d}/${s}/${s}_filtered_read_counts.txt" ]]          || return 1
  [[ -s "${d}/${s}/${s}_raw_coverage_per_smMIP.txt" ]]        || return 1
  [[ -s "${d}/${s}/${s}_UMI_usage_per_smMIP.txt" ]]           || return 1
  return 0
}

run_one() {
  local bam="$1"
  local base sample log
  base="$(basename "$bam")"
  sample="${base%%-*}"   # from {sample}-*-*--Sequences.coorSorted.bam
  if [[ -z "$sample" || "$sample" == "$base" ]]; then
    echo "⚠️  Could not derive sample from: $base" >&2
    return 0
  fi

  # Skip if outputs exist (unless FORCE_REDO=1)
  if sample_done "$sample" "$OUT_DIR" && [[ "$FORCE_REDO" != "1" ]]; then
    echo "⏭️  ${sample} already processed — skipping."
    return 0
  fi

  mkdir -p "${OUT_DIR}/${sample}"
  log="${LOG_DIR}/${sample}.smmip.log"

  echo "---------------------------------------------"
  echo "🔹 Processing: ${sample}"
  echo "    BAM:  ${bam}"
  echo "    OUT:  ${OUT_DIR}/${sample}"
  echo "    LOG:  ${log}"

  {
    echo "── $(date '+%F %T') :: starting ${sample}"
    echo "R_MAX_VSIZE=${R_MAX_VSIZE:-120GB} ; THREADS=${THREADS} ; MAPQ=${MAPQ} ; OVERLAP=${OVERLAP}"
  } | tee "$log"

  set +e
  R_MAX_VSIZE=120GB Rscript --no-restore --no-save \
    "$R_SCRIPT" \
    -b "$bam" \
    -p "$PANEL_FILE" \
    -s "$sample" \
    -o "$OUT_DIR" \
    -c "$CODE_DIR" \
    -t "$THREADS" \
    -O "$OVERLAP" \
    -M "$MAPQ"
  status=${PIPESTATUS[0]}
  set -e

  if [[ $status -ne 0 ]]; then
    echo "❌ R failed for ${sample} (exit $status). See: $log"
    return $status
  fi

  # Validate completion
  if ! sample_done "$sample" "$OUT_DIR"; then
    echo "❗ Expected outputs missing for ${sample}. See: $log"
    return 3
  fi

  echo "✅ Completed: ${sample}"
  # Hard reset R between samples to prevent memory creep
  clear_R_environment
  ensure_RMAX_in_Renviron
}

# ── Main ─────────────────────────────────────────────────────────────────────
ensure_RMAX_in_Renviron

# Find BAMs like *--Sequences.coorSorted.bam and run
while IFS= read -r -d '' bam; do
  run_one "$bam"
done < <(find "$BAM_DIR" -type f -name "*.filtered.bam" -print0)

echo "🎉 All eligible BAMs processed."
