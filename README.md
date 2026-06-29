# Food Allergy Variant Analysis Pipeline

Bioinformatics pipeline for identifying variants in STAT6 and FCER1A genes associated with food allergy.

## Workflow

1. Download SRA data
2. FASTQ Quality Control (FastQC)
3. Read trimming (fastp)
4. Sequence alignment (BWA)
5. BAM sorting and indexing
6. Variant calling (bcftools)
7. Variant extraction for STAT6 and FCER1A
8. SNP statistics and report generation

## Requirements

- SRA Toolkit
- FastQC
- fastp
- BWA
- samtools
- bcftools

## Usage

```bash
chmod +x pipeline.sh
./pipeline.sh
## Future Improvements

The following features will be incorporated in future versions:

- Variant annotation using ANNOVAR
- Variant annotation using Ensembl VEP
- Functional prediction using SIFT and PolyPhen
- Visualization of variant statistics
- Multi-sample comparative analysis
- Automated HTML report generation
- Workflow management using Snakemake
- Containerization using Docker
