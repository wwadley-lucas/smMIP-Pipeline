---
output:
  pdf_document: default
  html_document: default
---
# smMIP Analysis Pipeline Overview
## This pipeline utilizes the smMIPTools tool set that can be found at https://github.com/abelson-lab/smMIP-tools
## This pipeline impliments alignment, pre-filtering BAM files systems that may get bloated in R, parallel processing and sample/patient reports based on final smMIPTools output

## Quick Start

```bash
# Basic pipeline (Stages 1-6: FASTQ to mutation calls)
bash /path/to/smMIP-Pipeline/Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/experiment/RAW_FASTQ

# Full pipeline with CHIP analysis (Stages 1-8)
bash /path/to/smMIP-Pipeline/Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/experiment/RAW_FASTQ \
  --chip_analysis
```

## Pipeline Stages

| Stage | Script | Description |
|-------|--------|-------------|
| 1 | BWA.sh | Align paired-end reads to hg19 reference |
| 2 | filter_bam.sh | Filter for quality (MAPQ>=50, properly paired) |
| 3 | Read_Processing.sh | Map reads to smMIP probes, extract UMIs |
| 4 | Level_Bas_Calls.sh | Generate variant pileups at smMIP level |
| 5 | Annotate_SNVs.R | Annotate panel with gene/protein/COSMIC info |
| 6 | calling_mutations.R | Call mutations using statistical error modeling |
| 7 | CHIP_analysis.R | Population-level CHIP analysis (oncoplots, VAF heatmaps) |
| 8 | generate_patient_reports.R | Generate individual patient HTML/PDF reports |

## Command-Line Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--directory` | Base experiment directory (outputs go here) |
| `--fastq_dir` | Directory containing raw FASTQ files |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--panel` | `Universal/Myeloid_Panel_Targets.chr.txt` | Panel design file |
| `--ref` | `Universal/hg19/hg19.fa` | Reference genome |
| `--annotated_panel` | `Universal/annotated_Myeloid_Panel_Targets.txt` | Pre-annotated panel |
| `--config` | `{directory}/config.txt` | Sample configuration file |
| `--threads` | 6 | Threads per sample |
| `--parallel_jobs` | 4 | Samples to run in parallel |
| `--use_parallel` | false | Use parallel versions of scripts (flag) |
| `--force` | false | Force re-run all stages even if outputs exist |
| `--start` | 1 | Start at this stage (1-8) |
| `--stop` | 8 | Stop after this stage (1-8) |

### CHIP Analysis Parameters (Stages 7-8)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--chip_analysis` | - | Shorthand for --stop 8 (CHIP analysis stages) |
| `--chip_vaf` | 0.01 | VAF threshold for CHIP filtering |
| `--chip_pval` | 0.05 | P-value threshold for CHIP filtering |
| `--skip_reports` | false | Skip patient report generation (Stage 8) |
| `--skip_clinvar` | false | Skip ClinVar API queries (faster reports) |
| `--report_format` | html | Report format: html or pdf |

## Directory Structure

```
{--directory}/
├── bwa_out/                    # Stage 1 output
│   ├── *.bam
│   ├── *.bam.bai
│   ├── *.flagstat.txt
│   └── _logs/
├── filtered_bam/               # Stage 2 output
│   ├── *.filtered.bam
│   ├── *.filtered.bam.bai
│   └── _logs/
├── Read_Processing/            # Stage 3 output
│   ├── {sample}/
│   │   ├── {sample}_clean.bam
│   │   ├── {sample}_filtered_read_counts.txt
│   │   ├── {sample}_raw_coverage_per_smMIP.txt
│   │   └── {sample}_UMI_usage_per_smMIP.txt
│   ├── pileup/                 # Stage 4 output
│   │   ├── {sample}_raw_pileup.txt
│   │   └── {sample}_sscs_pileup.txt
│   └── _logs/
├── results/                    # Stage 6 output
│   └── called_mutations.txt
├── CHIP_analysis/              # Stages 7-8 output
│   ├── CHIP_oncoplot.pdf/png
│   ├── VAF_heatmap.pdf/png
│   ├── VAF_dotplot.pdf/png
│   ├── co_mutation_heatmap.pdf
│   ├── gene_mutation_summary.csv
│   ├── variant_annotation_summary.csv
│   ├── co_mutation_fisher_results.csv
│   └── patient_reports/
│       ├── {sample}_CHIP_Report.html
│       └── annotation_cache/
└── config.txt                  # Sample configuration
```

## Pipeline Flow Diagram

```
FASTQ files (*.R1*.fastq.gz, *.R2*.fastq.gz)
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 1: BWA Alignment                                      │
│ Script: BWA.sh / BWA_parallel.sh                            │
│ Input:  *R1*.fastq.gz, *R2*.fastq.gz                        │
│ Output: *.bam, *.bam.bai, *.flagstat.txt                    │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 2: BAM Filtering                                      │
│ Script: filter_bam.sh / filter_bam_parallel.sh              │
│ Input:  *.bam from Stage 1                                  │
│ Output: *.filtered.bam                                      │
│ Filter: MAPQ >= 50, properly paired (-f 2)                  │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 3: Read Processing / UMI Extraction                   │
│ Script: Read_Processing.sh / Read_Processing_parallel.sh    │
│ R Script: map_smMIPs_extract_UMIs.R                         │
│ Input:  *.filtered.bam + Panel_Targets.txt                  │
│ Output: Per sample directory containing:                    │
│         - {sample}_clean.bam                                │
│         - {sample}_filtered_read_counts.txt                 │
│         - {sample}_raw_coverage_per_smMIP.txt               │
│         - {sample}_UMI_usage_per_smMIP.txt                  │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 4: Pileup Generation (Level Base Calls)               │
│ Script: Level_Bas_Calls.sh / Level_Bas_Calls_parallel.sh    │
│ R Script: smMIP_level_raw_and_consensus_pileups.R           │
│ Input:  {sample}_clean.bam + Panel_Targets.txt              │
│ Output: - {sample}_raw_pileup.txt                           │
│         - {sample}_sscs_pileup.txt (consensus)              │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 5: Panel Annotation (ONE-TIME per panel)              │
│ R Script: Annotate_SNVs.R                                   │
│ Input:  Panel_Targets.txt                                   │
│ Output: annotated_Panel_Targets.txt                         │
│ Annotations: gene, protein, COSMIC, MAF, CADD scores        │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 6: Mutation Calling                                   │
│ R Script: calling_mutations.R                               │
│ Input:  - Pileup folder (*_raw_pileup.txt, *_sscs_pileup)   │
│         - Configuration file (sample info)                  │
│         - Annotated panel file                              │
│ Output: called_mutations.txt                                │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
  called_mutations.txt
    │
    ▼ (Optional: --chip_analysis)
┌─────────────────────────────────────────────────────────────┐
│ Stage 7: CHIP Population Analysis                           │
│ R Script: CHIP_analysis.R                                   │
│ Input:  called_mutations.txt                                │
│ Output: - CHIP_oncoplot.pdf/png (mutation landscape)        │
│         - VAF_heatmap.pdf/png (allele frequency matrix)     │
│         - VAF_dotplot.pdf/png (VAF distribution)            │
│         - co_mutation_heatmap.pdf (gene co-occurrence)      │
│         - gene_mutation_summary.csv                         │
│         - variant_annotation_summary.csv                    │
│         - co_mutation_fisher_results.csv                    │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│ Stage 8: Patient Report Generation                          │
│ R Script: generate_patient_reports.R                        │
│ Input:  called_mutations.txt                                │
│ Output: - patient_reports/{sample}_CHIP_Report.html         │
│         Per-patient clinical-style reports with:            │
│         - CHIP gene mutation summary                        │
│         - ClinVar annotations                               │
│         - Clinical significance assessment                  │
│         - VAF visualizations                                │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
  Patient Reports (HTML/PDF)
```

## Configuration File Format

The configuration file (`config.txt`) describes sample metadata for mutation calling.

### Format

Tab-delimited file with columns:

| Column | Description |
|--------|-------------|
| id | Sample ID (must match pileup filenames) |
| type | Sample type: `case` or `control` |
| replicate | Technical replicate pairing (or `NA`) |

### Example

```
id	type	replicate
Sample1	case	NA
Sample2	case	NA
Sample3	control	NA
Sample4	case	rep1
Sample5	case	rep1
```

### Auto-Generation

If `--config` is not provided, the pipeline will:
1. Scan the pileup directory for `*_raw_pileup.txt` files
2. Extract sample IDs from filenames
3. Generate a template `config.txt` with all samples as `case`
4. **Stop and prompt user to edit the file**
5. Re-run with `--start 6` to continue

## Mutation Calling Parameters

Key parameters for `calling_mutations.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-s` | Required | Path to pileup folder |
| `-f` | Required | Configuration file (sample layout) |
| `-a` | Required | Annotated panel file |
| `-b` | "sum" | Binomial error rate method ("sum" or "max") |
| `-g` | 10 | Min coverage for overlap (R1/R2) |
| `-m` | 0.001 | MAF cutoff to exclude known SNPs |
| `-v` | 0.05 | VAF cutoff for error modeling |
| `-p` | 0.05 | P-value cutoff (Bonferroni corrected) |
| `-d` | 2 | Fold-difference for cross-sample VAF |
| `-t` | 1 | Number of threads |

## Output: called_mutations.txt

Key columns in the final mutation calls:

| Column | Description |
|--------|-------------|
| sample_ID | Sample identifier |
| smMIP | smMIP probe name |
| chr | Chromosome |
| pos | Position (hg19) |
| ref | Reference allele |
| alt | Alternative allele |
| gene | Gene symbol |
| protein | Protein change (e.g., p.V600E) |
| cosmic | COSMIC ID if known |
| maf | Minor allele frequency (population) |
| variant_type | SNV, insertion, deletion |
| cadd_scaled | CADD pathogenicity score |
| P-value | Raw p-value |
| P-value.Bonferroni | Bonferroni-corrected p-value |
| non.ref.counts | Non-reference allele counts |
| total.depth | Total read depth |
| allele.frequency | Variant allele frequency |
| SSCS.* | Consensus-based metrics |
| flags | Quality flags |

## Output: CHIP Analysis (Stages 7-8)

### Stage 7 Outputs

| File | Description |
|------|-------------|
| `CHIP_oncoplot.pdf/png` | Mutation landscape visualization showing genes x samples |
| `VAF_heatmap.pdf/png` | Heatmap of variant allele frequencies |
| `VAF_dotplot.pdf/png` | Distribution of VAFs per gene |
| `co_mutation_heatmap.pdf` | Gene co-occurrence/mutual exclusivity analysis |
| `gene_mutation_summary.csv` | Frequency statistics per gene |
| `variant_annotation_summary.csv` | Per-variant details with protein changes |
| `co_mutation_fisher_results.csv` | Fisher's exact test results for co-mutation |

### Stage 8 Outputs

| File | Description |
|------|-------------|
| `patient_reports/{sample}_CHIP_Report.html` | Individual patient reports with: |
| | - Clinically significant mutations (VAF ≥ 1%) |
| | - Emerging clones (VAF < 1%) |
| | - CHIP gene clinical summaries |
| | - ClinVar pathogenicity annotations |
| | - VAF visualizations |

### CHIP Analysis Parameters

| Parameter | CLI Flag | Default | Description |
|-----------|----------|---------|-------------|
| VAF threshold | `--chip_vaf` | 0.01 | Minimum VAF for CHIP filtering |
| P-value threshold | `--chip_pval` | 0.05 | Maximum p-value for filtering |
| Skip ClinVar | `--skip_clinvar` | false | Use local knowledge base only |
| Report format | `--report_format` | html | Output format: html or pdf |

## Usage Examples

### Full Pipeline Run

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs
```

### Parallel Processing (Recommended for Large Datasets)

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs \
  --use_parallel \
  --parallel_jobs 4 \
  --threads 6
```

### Resume from Specific Stage

```bash
# After editing config.txt, run only mutation calling
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs \
  --start 6 --stop 6
```

### Force Re-run All Stages

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs \
  --force
```

### Run CHIP Analysis on Existing Mutations

```bash
# Run only CHIP stages (7-8) on existing called_mutations.txt
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --start 7 --stop 8
```

### Full Pipeline with CHIP Analysis

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs \
  --chip_analysis \
  --report_format html
```

### CHIP Analysis with Custom Thresholds

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --start 7 --stop 8 \
  --chip_vaf 0.02 \
  --chip_pval 0.01 \
  --skip_clinvar  # Faster, uses local knowledge base only
```

### Custom Panel and Reference

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/fastqs \
  --panel /path/to/custom_panel.txt \
  --ref /path/to/hg38.fa \
  --annotated_panel /path/to/annotated_panel.txt
```

## Resource Requirements

### Memory

| Stage | Memory per Sample | Notes |
|-------|-------------------|-------|
| 1 (BWA) | ~8 GB | Scales with reference genome size |
| 2 (Filter) | ~4 GB | I/O bound |
| 3 (ReadProc) | ~30 GB | R memory-intensive |
| 4 (Pileup) | ~25 GB | R memory-intensive |
| 5 (Annotate) | ~8 GB | Network-bound (API calls) |
| 6 (Mutations) | ~16 GB | Depends on cohort size |
| 7 (CHIP) | ~16 GB | Scales with cohort size |
| 8 (Reports) | ~8 GB | Per-sample, ClinVar API calls |

### Recommended Settings

| System RAM | --parallel_jobs | --threads |
|------------|-----------------|-----------|
| 32 GB | 2 | 4 |
| 64 GB | 3 | 6 |
| 128 GB | 4 | 6 |
| 256 GB | 8 | 6 |

## Dependencies

### Command-Line Tools

- bwa (or bwa-mem2)
- samtools
- GNU parallel (for --use_parallel)

### R Packages

Core Pipeline:
- optparse
- data.table
- parallel
- cellbaseR (for annotation)
- IRanges
- Rsamtools

CHIP Analysis (Stages 7-8):
- tidyverse
- ComplexHeatmap (Bioconductor)
- circlize
- RColorBrewer
- rmarkdown
- knitr
- kableExtra
- rentrez (for ClinVar API)

## Troubleshooting

### Out of Memory

- Reduce `--parallel_jobs`
- Ensure R_MAX_VSIZE is set (automatically handled by pipeline)

### Missing R1/R2 Pairs

- Check FASTQ naming convention matches `*R1*.fastq.gz` pattern
- Use `--r1_glob` and `--r2_token` to customize pattern

### Stage Skip Not Working

- Use `--force` to override skip logic
- Check that expected output files exist and are non-empty

### Config File Issues

- Ensure tab-delimited format
- Sample IDs must match pileup filenames exactly
- All samples need a type assignment (case or control)
