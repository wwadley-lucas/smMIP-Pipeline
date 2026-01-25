---
output:
  pdf_document: default
  html_document: default
---
# smMIP Analysis Pipeline Overview

## Quick Start

```bash
bash /Volumes/Seq_SSD/smMIP/Universal/Code/Master_pipeline.sh \
  --directory /Volumes/Seq_SSD/smMIP/KG001_01.22.25 \
  --fastq_dir /Volumes/Seq_SSD/smMIP/KG001_01.22.25/RAW_FASTQ
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
| `--start` | 1 | Start at this stage (1-6) |
| `--stop` | 6 | Stop after this stage (1-6) |

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

- optparse
- data.table
- parallel
- cellbaseR (for annotation)
- IRanges
- Rsamtools

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
