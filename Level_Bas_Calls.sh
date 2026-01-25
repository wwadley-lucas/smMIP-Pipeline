#!/usr/bin/env bash
set -euo pipefail

# ─── User paths ────────────────────────────────────────────────────────────────
READPROC_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/ReadProcessing"
PANEL_FILE="/Volumes/Seq_SSD/smMIP/10_23_25/Myeloid_Panel_Targets.chr.txt"
CODE_DIR="/Volumes/Seq_SSD/smMIP/10_23_25/code/smMIP-tools-main/R"
R_SCRIPT="${CODE_DIR}/smMIP_level_raw_and_consensus_pileups.R"

OUT_PILEUP="${READPROC_DIR}/pileup"
mkdir -p "${OUT_PILEUP}"

LOG_DIR="${OUT_PILEUP}/_logs"
mkdir -p "${LOG_DIR}"

# ─── Run settings (mirror your example) ───────────────────────────────────────
MMQ=20          # -m  minimum mapping quality
MBQ=25          # -q  minimum base quality
THREADS=4       # -t
# Optional pileup switches (keep defaults if you don't need them):
RANK="F"        # -r F/T  (report top alt only if 'T')
UMI="T"         # -u T/F  (include UMI columns)
FAM_SIZE=2      # -f
CONS_CUTOFF=0.7 # -v

# Rerun even if outputs exist: export FORCE_REDO=1
FORCE_REDO="${FORCE_REDO:-1}"

# ─── Helpers ──────────────────────────────────────────────────────────────────
ensure_RMAX_in_Renviron() {
  local renv="$HOME/.Renviron"
  touch "$renv"
  # remove any old setting, then set fresh
  if grep -q '^R_MAX_VSIZE=' "$renv"; then
    sed -i.bak '/^R_MAX_VSIZE=/d' "$renv"
  fi
  echo 'R_MAX_VSIZE=120GB' >> "$renv"
  export R_MAX_VSIZE="120GB"
  echo "ℹ️  R_MAX_VSIZE set to ${R_MAX_VSIZE}"
}

clear_R_environment() {
  # hard reset R between samples to avoid memory creep
  R_MAX_VSIZE=120GB Rscript --no-restore --no-save -e "rm(list=ls()); invisible(gc());"
}

# "Done" if we already have a non-empty raw or consensus pileup for the sample
sample_pileup_done() {
  local sample="$1"
  [[ -s "${OUT_PILEUP}/${sample}_raw_pileup.txt" ]] && return 0
  [[ -s "${OUT_PILEUP}/${sample}_consensus_pileup.txt" ]] && return 0
  return 1
}

run_one() {
  local sample_dir="$1"
  local sample bam log

  # sample name is the directory basename (e.g., ReadProcessing/10A/ → "10A")
  sample="$(basename "$sample_dir")"
  [[ -z "$sample" ]] && return 0

  # Prefer <SAMPLE>_clean.bam; otherwise use newest *clean.bam in the folder
  bam="${sample_dir%/}/${sample}_clean.bam"
  if [[ ! -f "$bam" ]]; then
    bam="$(find "$sample_dir" -type f -name "*clean.bam" -print0 \
           | xargs -0 ls -1t 2>/dev/null | head -n 1 || true)"
  fi
  if [[ -z "${bam:-}" || ! -f "$bam" ]]; then
    echo "⚠️  No clean BAM found for sample '${sample}' under: $sample_dir" >&2
    return 0
  fi

  # Skip if outputs already exist (unless FORCE_REDO=1)
  if sample_pileup_done "$sample" && [[ "$FORCE_REDO" != "1" ]]; then
    echo "⏭️  ${sample} already processed — skipping."
    return 0
  fi

  log="${LOG_DIR}/${sample}.pileup.log"
  echo "---------------------------------------------"
  echo "🔹 Processing: ${sample}"
  echo "    BAM: $bam"
  echo "    OUT: ${OUT_PILEUP}"
  echo "    LOG: ${log}"

  {
    echo "── $(date '+%F %T') :: starting ${sample}"
    echo "R_MAX_VSIZE=${R_MAX_VSIZE:-120GB} ; THREADS=${THREADS} ; MMQ=${MMQ} ; MBQ=${MBQ} ; RANK=${RANK} ; UMI=${UMI} ; FAM_SIZE=${FAM_SIZE} ; CONS_CUTOFF=${CONS_CUTOFF}"
  } | tee "$log"

  set +e
  R_MAX_VSIZE=120GB Rscript --no-restore --no-save "$R_SCRIPT" \
    -b "$bam" \
    -p "$PANEL_FILE" \
    -s "$sample" \
    -o "$OUT_PILEUP" \
    -c "$CODE_DIR" \
    -m "$MMQ" \
    -q "$MBQ" \
    -t "$THREADS" \
    -r "$RANK" \
    -u "$UMI" \
    -f "$FAM_SIZE" \
    -v "$CONS_CUTOFF" 2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}
  set -e

  if [[ $status -ne 0 ]]; then
    echo "❌ R failed for ${sample} (exit $status). See: $log"
    return $status
  fi

  # Validate completion
  if ! sample_pileup_done "$sample"; then
    echo "❗ Expected pileup outputs not found for ${sample}. See: $log"
    return 3
  fi

  echo "✅ Completed: ${sample}"
  clear_R_environment
  ensure_RMAX_in_Renviron
}

# ─── Main ─────────────────────────────────────────────────────────────────────
ensure_RMAX_in_Renviron

# Use find to be robust to odd names; only 1-level subdirs in READPROC_DIR
# (sample folders created by previous step).
while IFS= read -r -d '' d; do
  run_one "$d"
done < <(find "$READPROC_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

echo "🎉 All eligible samples processed."
