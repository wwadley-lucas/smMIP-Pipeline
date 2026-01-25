#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Master_pipeline.sh - Complete smMIP Analysis Pipeline
#
# Runs the entire smMIP analysis pipeline from raw FASTQ files to mutation calls
# =============================================================================

# ----------- Resource Path Resolution -----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # Universal/Code
UNIVERSAL_DIR="$(dirname "$SCRIPT_DIR")"                     # Universal

# ----------- Defaults -----------
DIRECTORY=""
FASTQ_DIR=""
BAM_DIR=""
PANEL="${UNIVERSAL_DIR}/Myeloid_Panel_Targets.chr.txt"
REF="${UNIVERSAL_DIR}/hg19/hg19.fa"
ANNOTATED_PANEL="${UNIVERSAL_DIR}/annotated_Myeloid_Panel_Targets.txt"
CONFIG=""
THREADS=6
PARALLEL_JOBS=4
USE_PARALLEL=false
FORCE=false
START_STAGE=1
STOP_STAGE=6

# R script directory
CODE_DIR="${UNIVERSAL_DIR}/Code/R"

# Pileup settings (can be customized)
MMQ=20
MBQ=25
OVERLAP=0.95
MAPQ=50

# Mutation calling settings
PVAL=0.05
VAF=0.05
MAF=0.001

# ----------- Usage -----------
usage() {
  cat <<EOF
Usage:
  $0 --directory DIR --fastq_dir DIR [options]
  $0 --directory DIR --bam_dir DIR --start 3 [options]

Required:
  --directory DIR        Base experiment directory (outputs go here)
  --fastq_dir DIR        Directory containing raw FASTQ files (required for stages 1-2)

Optional:
  --bam_dir DIR          Directory with pre-aligned BAM files (use with --start 3)
  --panel FILE           Panel design file (default: Universal/Myeloid_Panel_Targets.chr.txt)
  --ref FILE             Reference genome (default: Universal/hg19/hg19.fa)
  --annotated_panel FILE Pre-annotated panel (default: Universal/annotated_Myeloid_Panel_Targets.txt)
  --config FILE          Sample configuration file (default: {directory}/config.txt)
  --threads N            Threads per sample (default: 6)
  --parallel_jobs N      Samples to run in parallel (default: 4)
  --use_parallel         Use parallel versions of scripts (flag)
  --force                Force re-run all stages even if outputs exist
  --start N              Start at this stage (1-6, default: 1)
  --stop N               Stop after this stage (1-6, default: 6)
  -h, --help             Show this help message

Stages:
  1 - BWA alignment
  2 - BAM filtering
  3 - Read processing / UMI extraction
  4 - Pileup generation (Level base calls)
  5 - Panel annotation (skipped if annotated panel exists)
  6 - Mutation calling

Example:
  $0 --directory /Volumes/Seq_SSD/smMIP/KG001_01.22.25 \\
     --fastq_dir /Volumes/Seq_SSD/smMIP/KG001_01.22.25/RAW_FASTQ

Output Structure:
  {directory}/
  ├── bwa_out/           # Stage 1: Aligned BAMs
  ├── filtered_bam/      # Stage 2: Filtered BAMs
  ├── Read_Processing/   # Stage 3: Processed reads
  │   └── pileup/        # Stage 4: Pileup files
  └── results/           # Stage 6: Mutation calls
EOF
}

# ----------- Helper Functions -----------
say() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
error() { say "ERROR: $*" >&2; exit 1; }
warn() { say "WARNING: $*" >&2; }

check_file() {
  [[ -f "$1" ]] || error "File not found: $1"
}

check_dir() {
  [[ -d "$1" ]] || error "Directory not found: $1"
}

# Count files matching a pattern
count_files() {
  local dir="$1"
  local pattern="$2"
  find "$dir" -maxdepth 1 -name "$pattern" -type f 2>/dev/null | wc -l | tr -d ' '
}

# Count FASTQ samples (R1 files)
count_fastq_samples() {
  count_files "$FASTQ_DIR" "*R1*.fastq.gz"
}

# Check if stage should be skipped
should_skip_stage() {
  local stage="$1"
  [[ "$FORCE" == true ]] && return 1  # Never skip if force is enabled

  case "$stage" in
    1)  # BWA - skip if BAM count matches FASTQ count
        local expected=$(count_fastq_samples)
        local actual=$(count_files "$BWA_OUT" "*.bam")
        [[ "$actual" -ge "$expected" && "$expected" -gt 0 ]]
        ;;
    2)  # Filter - skip if filtered BAM count matches BWA output
        local expected=$(count_files "$BWA_OUT" "*.bam")
        local actual=$(count_files "$FILTERED_BAM" "*.filtered.bam")
        [[ "$actual" -ge "$expected" && "$expected" -gt 0 ]]
        ;;
    3)  # Read Processing - skip if clean BAMs exist for all samples
        # Count input BAMs (handle both .filtered.bam and .bam patterns)
        local expected=$(count_files "$FILTERED_BAM" "*.filtered.bam")
        [[ "$expected" -eq 0 ]] && expected=$(count_files "$FILTERED_BAM" "*.bam")
        local actual=0
        for d in "$READ_PROC_DIR"/*/; do
          [[ -d "$d" ]] || continue
          local sample=$(basename "$d")
          [[ -f "${d}${sample}_clean.bam" ]] && ((actual++)) || true
        done
        [[ "$actual" -ge "$expected" && "$expected" -gt 0 ]]
        ;;
    4)  # Pileup - skip if pileup files exist for all processed samples
        local actual=0
        for d in "$READ_PROC_DIR"/*/; do
          [[ -d "$d" ]] || continue
          local sample=$(basename "$d")
          [[ "$sample" == "pileup" || "$sample" == "_logs" ]] && continue
          [[ -f "$PILEUP_DIR/${sample}_raw_pileup.txt" ]] && ((actual++)) || true
        done
        local expected=0
        for d in "$READ_PROC_DIR"/*/; do
          [[ -d "$d" ]] || continue
          local sample=$(basename "$d")
          [[ "$sample" == "pileup" || "$sample" == "_logs" ]] && continue
          ((expected++)) || true
        done
        [[ "$actual" -ge "$expected" && "$expected" -gt 0 ]]
        ;;
    5)  # Annotate - skip if annotated panel exists
        [[ -f "$ANNOTATED_PANEL" ]]
        ;;
    6)  # Mutations - never skip (always re-run)
        return 1
        ;;
    *)
        return 1
        ;;
  esac
}

# Generate config file template
generate_config() {
  local config_file="$1"
  local pileup_dir="$2"

  say "Generating config file template: $config_file"

  # Header
  echo -e "id\ttype\treplicate" > "$config_file"

  # Extract sample IDs from pileup files
  for f in "$pileup_dir"/*_raw_pileup.txt; do
    [[ -f "$f" ]] || continue
    local sample=$(basename "$f" _raw_pileup.txt)
    echo -e "${sample}\tcase\tNA" >> "$config_file"
  done

  say "Config template generated with $(( $(wc -l < "$config_file") - 1 )) samples"
  warn "Please review and edit $config_file to set sample types (case/control) before running Stage 6"
}

# ----------- Parse Arguments -----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --directory) DIRECTORY="$2"; shift 2;;
    --fastq_dir) FASTQ_DIR="$2"; shift 2;;
    --bam_dir) BAM_DIR="$2"; shift 2;;
    --panel) PANEL="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --annotated_panel) ANNOTATED_PANEL="$2"; shift 2;;
    --config) CONFIG="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --parallel_jobs) PARALLEL_JOBS="$2"; shift 2;;
    --use_parallel) USE_PARALLEL=true; shift;;
    --force) FORCE=true; shift;;
    --start) START_STAGE="$2"; shift 2;;
    --stop) STOP_STAGE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) error "Unknown argument: $1";;
  esac
done

# ----------- Validate Required Arguments -----------
[[ -z "$DIRECTORY" ]] && { usage; error "--directory is required"; }

# Validate stage range
[[ "$START_STAGE" -lt 1 || "$START_STAGE" -gt 6 ]] && error "--start must be between 1 and 6"
[[ "$STOP_STAGE" -lt 1 || "$STOP_STAGE" -gt 6 ]] && error "--stop must be between 1 and 6"
[[ "$START_STAGE" -gt "$STOP_STAGE" ]] && error "--start cannot be greater than --stop"

# Validate inputs based on start stage
if [[ "$START_STAGE" -lt 3 ]]; then
  [[ -z "$FASTQ_DIR" ]] && { usage; error "--fastq_dir is required for stages 1-2"; }
  check_dir "$FASTQ_DIR"
else
  # Starting from stage 3+ - need either BAM_DIR or we'll check for existing filtered_bam later
  if [[ -n "$BAM_DIR" ]]; then
    check_dir "$BAM_DIR"
  elif [[ -z "$FASTQ_DIR" ]]; then
    # Neither BAM_DIR nor FASTQ_DIR provided - check if filtered_bam exists
    if [[ ! -d "${DIRECTORY}/filtered_bam" ]]; then
      usage
      error "Starting at stage 3 requires --bam_dir or existing filtered_bam/ directory"
    fi
  fi
fi

# ----------- Validate Paths -----------
check_file "$PANEL"
check_file "$REF"
[[ -d "$CODE_DIR" ]] || error "R code directory not found: $CODE_DIR"

# ----------- Derive Output Directories -----------
BWA_OUT="${DIRECTORY}/bwa_out"
FILTERED_BAM="${DIRECTORY}/filtered_bam"
READ_PROC_DIR="${DIRECTORY}/Read_Processing"
PILEUP_DIR="${READ_PROC_DIR}/pileup"
RESULTS_DIR="${DIRECTORY}/results"

# Use external BAM_DIR as input for Stage 3 if provided
if [[ -n "$BAM_DIR" ]]; then
  FILTERED_BAM="$BAM_DIR"
fi

# Set default config if not provided
[[ -z "$CONFIG" ]] && CONFIG="${DIRECTORY}/config.txt"

# ----------- Create Directories -----------
mkdir -p "$BWA_OUT" "$FILTERED_BAM" "$READ_PROC_DIR" "$PILEUP_DIR" "$RESULTS_DIR"

# ----------- Print Configuration -----------
say "============================================="
say "       smMIP Analysis Pipeline"
say "============================================="
say "Directory:        $DIRECTORY"
[[ -n "$FASTQ_DIR" ]] && say "FASTQ Dir:        $FASTQ_DIR"
[[ -n "$BAM_DIR" ]] && say "BAM Dir:          $BAM_DIR"
say "Panel:            $PANEL"
say "Reference:        $REF"
say "Annotated Panel:  $ANNOTATED_PANEL"
say "Config:           $CONFIG"
say "Threads:          $THREADS"
say "Parallel Jobs:    $PARALLEL_JOBS"
say "Use Parallel:     $USE_PARALLEL"
say "Force Rerun:      $FORCE"
say "Stages:           $START_STAGE to $STOP_STAGE"
say "============================================="

# =============================================================================
# STAGE 1: BWA Alignment
# =============================================================================
if [[ "$START_STAGE" -le 1 && "$STOP_STAGE" -ge 1 ]]; then
  say ""
  say "========== STAGE 1: BWA Alignment =========="

  if should_skip_stage 1; then
    say "Skipping Stage 1 (outputs exist). Use --force to re-run."
  else
    if [[ "$USE_PARALLEL" == true ]]; then
      say "Running BWA_parallel.sh..."
      bash "${SCRIPT_DIR}/BWA_parallel.sh" \
        --fastq_dir "$FASTQ_DIR" \
        --ref "$REF" \
        --out "$BWA_OUT" \
        --threads "$THREADS" \
        --parallel "$PARALLEL_JOBS"
    else
      say "Running BWA.sh..."
      bash "${SCRIPT_DIR}/BWA.sh" \
        --fastq_dir "$FASTQ_DIR" \
        --ref "$REF" \
        --out "$BWA_OUT" \
        --threads "$THREADS"
    fi
    say "Stage 1 complete."
  fi
fi

# =============================================================================
# STAGE 2: BAM Filtering
# =============================================================================
if [[ "$START_STAGE" -le 2 && "$STOP_STAGE" -ge 2 ]]; then
  say ""
  say "========== STAGE 2: BAM Filtering =========="

  if should_skip_stage 2; then
    say "Skipping Stage 2 (outputs exist). Use --force to re-run."
  else
    if [[ "$USE_PARALLEL" == true ]]; then
      say "Running filter_bam_parallel.sh..."
      bash "${SCRIPT_DIR}/filter_bam_parallel.sh" \
        --bam_dir "$BWA_OUT" \
        --out "$FILTERED_BAM" \
        --threads "$THREADS" \
        --parallel "$PARALLEL_JOBS" \
        --mapq "$MAPQ"
    else
      say "Running filter_bam.sh..."
      bash "${SCRIPT_DIR}/filter_bam.sh" \
        --bam_dir "$BWA_OUT" \
        --out "$FILTERED_BAM" \
        --threads "$THREADS" \
        --mapq "$MAPQ"
    fi
    say "Stage 2 complete."
  fi
fi

# =============================================================================
# STAGE 3: Read Processing / UMI Extraction
# =============================================================================
if [[ "$START_STAGE" -le 3 && "$STOP_STAGE" -ge 3 ]]; then
  say ""
  say "========== STAGE 3: Read Processing =========="

  if should_skip_stage 3; then
    say "Skipping Stage 3 (outputs exist). Use --force to re-run."
  else
    if [[ "$USE_PARALLEL" == true ]]; then
      say "Running Read_Processing_parallel.sh..."
      # Set FORCE_REDO based on --force flag
      FORCE_FLAG=0
      [[ "$FORCE" == true ]] && FORCE_FLAG=1

      FORCE_REDO=$FORCE_FLAG bash "${SCRIPT_DIR}/Read_Processing_parallel.sh" \
        --bam_dir "$FILTERED_BAM" \
        --panel "$PANEL" \
        --code_dir "$CODE_DIR" \
        --out "$READ_PROC_DIR" \
        --threads "$THREADS" \
        --parallel "$PARALLEL_JOBS" \
        --overlap "$OVERLAP" \
        --mapq "$MAPQ"
    else
      say "Running Read_Processing.sh (sequential)..."
      # The sequential script uses hardcoded paths, so we need to set env vars
      # and run with modified variables
      FORCE_FLAG=0
      [[ "$FORCE" == true ]] && FORCE_FLAG=1

      # Create a temporary wrapper to run with correct paths
      BAM_DIR="$FILTERED_BAM" \
      PANEL_FILE="$PANEL" \
      CODE_DIR="$CODE_DIR" \
      OUT_DIR="$READ_PROC_DIR" \
      THREADS="$THREADS" \
      OVERLAP="$OVERLAP" \
      MAPQ="$MAPQ" \
      FORCE_REDO="$FORCE_FLAG" \
      R_SCRIPT="${CODE_DIR}/map_smMIPs_extract_UMIs.R" \
      bash -c '
        set -euo pipefail

        mkdir -p "${OUT_DIR}/_logs"
        LOG_DIR="${OUT_DIR}/_logs"

        sample_done() {
          local s="$1" d="$2"
          [[ -s "${d}/${s}/${s}_clean.bam" ]] || return 1
          [[ -s "${d}/${s}/${s}_filtered_read_counts.txt" ]] || return 1
          return 0
        }

        # Find BAM files - try .filtered.bam first, then fall back to .bam
        shopt -s nullglob
        bam_files=("$BAM_DIR"/*.filtered.bam)
        bam_suffix=".filtered.bam"
        if [[ ${#bam_files[@]} -eq 0 ]]; then
          bam_files=("$BAM_DIR"/*.bam)
          bam_suffix=".bam"
        fi

        for bam in "${bam_files[@]}"; do
          base="$(basename "$bam")"
          sample="${base%$bam_suffix}"

          if sample_done "$sample" "$OUT_DIR" && [[ "$FORCE_REDO" != "1" ]]; then
            echo "[$(date "+%Y-%m-%d %H:%M:%S")] $sample already processed - skipping."
            continue
          fi

          mkdir -p "${OUT_DIR}/${sample}"
          log="${LOG_DIR}/${sample}.smmip.log"

          echo "[$(date "+%Y-%m-%d %H:%M:%S")] Processing: ${sample}"

          R_MAX_VSIZE=120GB Rscript --no-restore --no-save \
            "$R_SCRIPT" \
            -b "$bam" \
            -p "$PANEL_FILE" \
            -s "$sample" \
            -o "$OUT_DIR" \
            -c "$CODE_DIR" \
            -t "$THREADS" \
            -O "$OVERLAP" \
            -M "$MAPQ" 2>&1 | tee "$log"

          echo "[$(date "+%Y-%m-%d %H:%M:%S")] Completed: ${sample}"
        done
      '
    fi
    say "Stage 3 complete."
  fi
fi

# =============================================================================
# STAGE 4: Pileup Generation (Level Base Calls)
# =============================================================================
if [[ "$START_STAGE" -le 4 && "$STOP_STAGE" -ge 4 ]]; then
  say ""
  say "========== STAGE 4: Pileup Generation =========="

  if should_skip_stage 4; then
    say "Skipping Stage 4 (outputs exist). Use --force to re-run."
  else
    if [[ "$USE_PARALLEL" == true ]]; then
      say "Running Level_Bas_Calls_parallel.sh..."
      FORCE_FLAG=0
      [[ "$FORCE" == true ]] && FORCE_FLAG=1

      FORCE_REDO=$FORCE_FLAG bash "${SCRIPT_DIR}/Level_Bas_Calls_parallel.sh" \
        --readproc_dir "$READ_PROC_DIR" \
        --panel "$PANEL" \
        --code_dir "$CODE_DIR" \
        --out "$PILEUP_DIR" \
        --threads "$THREADS" \
        --parallel "$PARALLEL_JOBS" \
        --mmq "$MMQ" \
        --mbq "$MBQ"
    else
      say "Running Level_Bas_Calls.sh (sequential)..."
      FORCE_FLAG=0
      [[ "$FORCE" == true ]] && FORCE_FLAG=1

      READPROC_DIR="$READ_PROC_DIR" \
      PANEL_FILE="$PANEL" \
      CODE_DIR="$CODE_DIR" \
      OUT_PILEUP="$PILEUP_DIR" \
      THREADS="$THREADS" \
      MMQ="$MMQ" \
      MBQ="$MBQ" \
      FORCE_REDO="$FORCE_FLAG" \
      R_SCRIPT="${CODE_DIR}/smMIP_level_raw_and_consensus_pileups.R" \
      bash -c '
        set -euo pipefail

        mkdir -p "${OUT_PILEUP}/_logs"
        LOG_DIR="${OUT_PILEUP}/_logs"

        RANK="F"
        UMI="T"
        FAM_SIZE=2
        CONS_CUTOFF=0.7

        sample_pileup_done() {
          local sample="$1"
          [[ -s "${OUT_PILEUP}/${sample}_raw_pileup.txt" ]] && return 0
          [[ -s "${OUT_PILEUP}/${sample}_consensus_pileup.txt" ]] && return 0
          return 1
        }

        for sample_dir in "$READPROC_DIR"/*/; do
          [[ -d "$sample_dir" ]] || continue
          sample="$(basename "$sample_dir")"
          [[ "$sample" == "pileup" || "$sample" == "_logs" ]] && continue

          bam="${sample_dir}${sample}_clean.bam"
          [[ -f "$bam" ]] || continue

          if sample_pileup_done "$sample" && [[ "$FORCE_REDO" != "1" ]]; then
            echo "[$(date "+%Y-%m-%d %H:%M:%S")] $sample already processed - skipping."
            continue
          fi

          log="${LOG_DIR}/${sample}.pileup.log"
          echo "[$(date "+%Y-%m-%d %H:%M:%S")] Processing: ${sample}"

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
            -v "$CONS_CUTOFF" 2>&1 | tee "$log"

          echo "[$(date "+%Y-%m-%d %H:%M:%S")] Completed: ${sample}"
        done
      '
    fi
    say "Stage 4 complete."
  fi
fi

# =============================================================================
# STAGE 5: Panel Annotation
# =============================================================================
if [[ "$START_STAGE" -le 5 && "$STOP_STAGE" -ge 5 ]]; then
  say ""
  say "========== STAGE 5: Panel Annotation =========="

  if should_skip_stage 5; then
    say "Skipping Stage 5 (annotated panel exists): $ANNOTATED_PANEL"
  else
    say "Running Annotate_SNVs.R..."
    say "This may take a while for first-time annotation..."

    # Determine output location - same directory as panel file
    ANNOTATION_OUTPUT_DIR="$(dirname "$PANEL")"

    R_MAX_VSIZE=120GB Rscript --no-restore --no-save \
      "${CODE_DIR}/Annotate_SNVs.R" \
      -p "$PANEL" \
      -c "$CODE_DIR" \
      -t "$THREADS" \
      -g "GRCh37"

    # The R script outputs to the panel file's directory with "annotated_" prefix
    EXPECTED_ANNOTATED="${ANNOTATION_OUTPUT_DIR}/annotated_$(basename "$PANEL")"
    if [[ -f "$EXPECTED_ANNOTATED" ]]; then
      ANNOTATED_PANEL="$EXPECTED_ANNOTATED"
      say "Annotation complete: $ANNOTATED_PANEL"
    else
      warn "Annotation output not found at expected location: $EXPECTED_ANNOTATED"
    fi

    say "Stage 5 complete."
  fi
fi

# =============================================================================
# STAGE 6: Mutation Calling
# =============================================================================
if [[ "$START_STAGE" -le 6 && "$STOP_STAGE" -ge 6 ]]; then
  say ""
  say "========== STAGE 6: Mutation Calling =========="

  # Check if config exists, generate template if not
  if [[ ! -f "$CONFIG" ]]; then
    generate_config "$CONFIG" "$PILEUP_DIR"
    say ""
    say "============================================="
    say "ACTION REQUIRED: Config file generated"
    say "============================================="
    say "Please edit: $CONFIG"
    say "Set the 'type' column to 'case' or 'control' for each sample"
    say "Then re-run the pipeline with --start 6"
    say "============================================="
    exit 0
  fi

  # Validate required files
  check_file "$CONFIG"
  check_file "$ANNOTATED_PANEL"

  say "Running calling_mutations.R..."

  R_MAX_VSIZE=120GB Rscript --no-restore --no-save \
    "${CODE_DIR}/calling_mutations.R" \
    -s "$PILEUP_DIR" \
    -f "$CONFIG" \
    -a "$ANNOTATED_PANEL" \
    -c "$CODE_DIR" \
    -t "$THREADS" \
    -p "$PVAL" \
    -v "$VAF" \
    -m "$MAF" \
    -o "$RESULTS_DIR"

  say "Stage 6 complete."

  # Check for output
  if [[ -f "${RESULTS_DIR}/called_mutations.txt" ]]; then
    MUTATION_COUNT=$(( $(wc -l < "${RESULTS_DIR}/called_mutations.txt") - 1 ))
    say "Output: ${RESULTS_DIR}/called_mutations.txt ($MUTATION_COUNT mutations)"
  fi
fi

# =============================================================================
# Pipeline Complete
# =============================================================================
say ""
say "============================================="
say "       Pipeline Complete!"
say "============================================="
say "Results directory: $RESULTS_DIR"

# List outputs
say ""
say "Output files:"
[[ -d "$BWA_OUT" ]] && say "  - BWA alignments:    $BWA_OUT ($(count_files "$BWA_OUT" "*.bam") BAMs)"
[[ -d "$FILTERED_BAM" ]] && say "  - Filtered BAMs:     $FILTERED_BAM ($(count_files "$FILTERED_BAM" "*.filtered.bam") BAMs)"
[[ -d "$PILEUP_DIR" ]] && say "  - Pileup files:      $PILEUP_DIR ($(count_files "$PILEUP_DIR" "*_raw_pileup.txt") samples)"
[[ -f "${RESULTS_DIR}/called_mutations.txt" ]] && say "  - Mutations:         ${RESULTS_DIR}/called_mutations.txt"

say ""
say "Done."
