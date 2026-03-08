# smMIP Targeted Sequencing Pipeline

An end-to-end analysis pipeline for single-molecule Molecular Inversion Probe (smMIP) targeted sequencing data. Takes paired-end FASTQ files through alignment, UMI-aware read processing, variant pileup generation, statistical mutation calling, and optional CHIP (Clonal Hematopoiesis of Indeterminate Potential) population analysis with per-patient clinical reports.

Built on top of [smMIP-tools](https://github.com/abelson-lab/smMIP-tools) (Abelson Lab, OICR) with added stages for alignment, BAM pre-filtering, parallel processing, and downstream CHIP analysis.

## Prerequisites

### Command-Line Tools

| Tool | Version | Purpose |
|------|---------|---------|
| BWA (or bwa-mem2) | 0.7.17+ | Read alignment to reference genome |
| samtools | 1.10+ | BAM filtering, sorting, indexing |
| GNU parallel | (optional) | Parallel sample processing (`--use_parallel`) |
| R | 4.4+ | All R-based analysis stages |

### R Packages

**Core pipeline (Stages 3-6):**

- optparse
- data.table
- parallel
- Rsamtools (Bioconductor)
- IRanges (Bioconductor)
- cellbaseR (Bioconductor) -- used by Stage 5 annotation

**CHIP analysis (Stages 7-8):**

- tidyverse
- ComplexHeatmap (Bioconductor)
- circlize
- RColorBrewer
- rmarkdown, knitr, kableExtra -- patient report rendering
- rentrez -- ClinVar API queries

Run `Rscript R/CHIP/install_deps.R` to install CHIP-stage dependencies automatically.

### Reference Data

- hg19 reference genome (indexed for BWA)
- smMIP panel design file (MIPgen format, tab-delimited)

## Pipeline Stages

| Stage | Name | Script(s) | Description |
|-------|------|-----------|-------------|
| 1 | BWA Alignment | `BWA.sh` / `BWA_parallel.sh` | Align paired-end FASTQ reads to hg19 with BWA-MEM |
| 2 | BAM Filtering | `filter_bam.sh` / `filter_bam_parallel.sh` | Filter for MAPQ >= 50 and properly paired reads |
| 3 | Read Processing | `Read_Processing.sh` / `Read_Processing_parallel.sh` | Map reads to smMIP probes and extract UMIs (`map_smMIPs_extract_UMIs.R`) |
| 4 | Pileup Generation | `Level_Bas_Calls.sh` / `Level_Bas_Calls_parallel.sh` | Generate raw and consensus (SSCS) pileups per smMIP (`smMIP_level_raw_and_consensus_pileups.R`) |
| 5 | Panel Annotation | `Annotate_SNVs.R` | Annotate panel positions with gene, protein, COSMIC, MAF, and CADD scores (one-time per panel) |
| 6 | Mutation Calling | `calling_mutations.R` | Call mutations using binomial error modeling with Bonferroni correction |
| 7 | CHIP Analysis | `CHIP_analysis.R` | Population-level oncoplots, VAF heatmaps, co-mutation analysis |
| 8 | Patient Reports | `generate_patient_reports.R` | Per-patient HTML/PDF reports with ClinVar annotations |

## Usage

### Basic Run (Stages 1-6)

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/experiment/RAW_FASTQ
```

### Full Pipeline with CHIP Analysis (Stages 1-8)

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/experiment/RAW_FASTQ \
  --chip_analysis
```

### Parallel Processing (recommended for large cohorts)

```bash
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --fastq_dir /path/to/experiment/RAW_FASTQ \
  --use_parallel --parallel_jobs 4 --threads 6
```

### Resume from a Specific Stage

```bash
# Re-run only mutation calling after editing config.txt
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --start 6 --stop 6

# Run CHIP stages on existing called_mutations.txt
bash Master_pipeline.sh \
  --directory /path/to/experiment \
  --start 7 --stop 8
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

### Key Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--directory` | (required) | Base experiment directory for all outputs |
| `--fastq_dir` | (required for stages 1-2) | Directory containing `*R1*.fastq.gz` / `*R2*.fastq.gz` |
| `--bam_dir` | -- | Pre-aligned BAMs (use with `--start 3`) |
| `--panel` | `Universal/Myeloid_Panel_Targets.chr.txt` | smMIP panel design file |
| `--ref` | `Universal/hg19/hg19.fa` | BWA-indexed reference genome |
| `--config` | `{directory}/config.txt` | Sample configuration (case/control assignments) |
| `--threads` | 2 | Threads per sample |
| `--parallel_jobs` | 4 | Samples to process in parallel |
| `--use_parallel` | false | Enable parallel versions of shell scripts |
| `--force` | false | Re-run stages even if outputs already exist |
| `--start` / `--stop` | 1 / 8 | Stage range to execute |
| `--chip_analysis` | -- | Shorthand for `--stop 8` |
| `--chip_vaf` | 0.01 | VAF threshold for CHIP filtering |
| `--chip_pval` | 0.05 | P-value threshold for CHIP filtering |
| `--skip_reports` | false | Skip Stage 8 patient report generation |
| `--skip_clinvar` | false | Skip ClinVar API queries (use local knowledge base only) |
| `--report_format` | html | Patient report format: `html` or `pdf` |

## Input Format

### FASTQ Files

Paired-end FASTQ files with naming convention `*R1*.fastq.gz` and `*R2*.fastq.gz`. Place all samples in a single directory passed via `--fastq_dir`.

### Panel Design File

Tab-delimited MIPgen-format file describing smMIP probe coordinates and sequences. Must contain columns: `mip_name`, `chr`, `ext_probe_start`, `ext_probe_stop`, `lig_probe_start`, `lig_probe_stop`, `ext_probe_sequence`, `lig_probe_sequence`, `mip_sequence`, `scan_target_sequence`, `mip_scan_start_position`, `mip_scan_stop_position`, `probe_strand`.

### Sample Configuration (`config.txt`)

Tab-delimited file with three columns:

```
id	type	replicate
Sample1	case	NA
Sample2	case	NA
Sample3	control	NA
Sample4	case	rep1
Sample5	case	rep1
```

- **id**: Must match pileup filenames exactly (derived from FASTQ sample names).
- **type**: `case` or `control`. Controls are used for error modeling.
- **replicate**: Technical replicate grouping, or `NA`.

If `--config` is not provided and no `config.txt` exists, the pipeline auto-generates a template at Stage 6, populates it with all samples as `case`, and pauses for the user to review before continuing.

## Output

```
{--directory}/
├── bwa_out/                    # Stage 1: aligned BAMs + flagstat
├── filtered_bam/               # Stage 2: quality-filtered BAMs
├── Read_Processing/            # Stage 3: per-sample UMI-processed reads
│   ├── {sample}/
│   │   ├── {sample}_clean.bam
│   │   ├── {sample}_filtered_read_counts.txt
│   │   ├── {sample}_raw_coverage_per_smMIP.txt
│   │   └── {sample}_UMI_usage_per_smMIP.txt
│   └── pileup/                 # Stage 4: variant pileups
│       ├── {sample}_raw_pileup.txt
│       └── {sample}_sscs_pileup.txt
├── results/                    # Stage 6: final mutation calls
│   └── called_mutations.txt
├── CHIP_analysis/              # Stages 7-8
│   ├── CHIP_oncoplot.pdf/png
│   ├── VAF_heatmap.pdf/png
│   ├── VAF_dotplot.pdf/png
│   ├── co_mutation_heatmap.pdf
│   ├── gene_mutation_summary.csv
│   ├── variant_annotation_summary.csv
│   ├── co_mutation_fisher_results.csv
│   └── patient_reports/
│       └── {sample}_CHIP_Report.html
└── config.txt
```

### `called_mutations.txt` Columns

| Column | Description |
|--------|-------------|
| sample_ID | Sample identifier |
| chr, pos, ref, alt | Genomic coordinates and alleles (hg19) |
| gene | Gene symbol |
| protein | Protein change (e.g., p.V600E) |
| cosmic | COSMIC ID if annotated |
| maf | Population minor allele frequency |
| variant_type | SNV, insertion, or deletion |
| cadd_scaled | CADD pathogenicity score |
| P-value / P-value.Bonferroni | Raw and corrected p-values |
| non.ref.counts / total.depth | Variant and total read counts |
| allele.frequency | Variant allele frequency |
| SSCS.* | Consensus sequence-based metrics |
| flags | Quality warning flags (empty = passed) |

## Known Limitations

- **Reference genome**: The pipeline defaults to hg19. Using hg38 requires supplying `--ref` and `--panel` files built for GRCh38, and the annotation stage (Stage 5) must be re-run with `-g GRCh38`.
- **Memory usage**: Stages 3 and 4 are R memory-intensive (~25-30 GB per sample). Reduce `--parallel_jobs` on systems with limited RAM.
- **Panel annotation (Stage 5)**: Uses the CellBase API for variant annotation, which requires network access and can be slow on first run. Results are cached via the annotated panel file.
- **ClinVar queries (Stage 8)**: Patient reports query the NCBI ClinVar API in real time. Use `--skip_clinvar` for faster report generation with local knowledge base only.
- **FASTQ naming**: Expects `*R1*.fastq.gz` and `*R2*.fastq.gz` naming. Non-standard naming will cause Stage 1 to find no input files.
- **Config auto-generation**: The auto-generated `config.txt` marks all samples as `case`. Users must manually set control samples before running Stage 6.
- **Parallel scripts**: The `--use_parallel` flag requires GNU `parallel` to be installed and on the PATH.
- **Single-panel design**: The pipeline processes one panel design per run. Multi-panel experiments require separate pipeline invocations.

## License

MIT License. Copyright (c) 2021 Ontario Institute for Cancer Research. See [LICENSE](LICENSE) for details.
