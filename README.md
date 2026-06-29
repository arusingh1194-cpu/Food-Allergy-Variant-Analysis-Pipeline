# Food Allergy Variant Analysis Pipeline

## Overview

This project presents an automated bioinformatics pipeline for identifying and analyzing genetic variants associated with food allergy susceptibility, focusing on the **STAT6** and **FCER1A** genes.

The workflow performs:

* SRA data download
* Quality control
* Read trimming
* Sequence alignment
* Variant calling
* Variant annotation
* Functional prediction
* Gene-wise variant extraction
* Visualization
* Automated report generation

---

## Workflow

```text
SRA Download
      ↓
FASTQ Quality Control (FastQC)
      ↓
Read Trimming (fastp)
      ↓
Post-trimming Quality Check (FastQC)
      ↓
Reference Genome Indexing (BWA)
      ↓
Read Alignment (BWA-MEM)
      ↓
BAM Sorting and Indexing (SAMtools)
      ↓
Variant Calling (BCFtools)
      ↓
Variant Annotation (Ensembl VEP)
      ↓
Functional Prediction
(SIFT + PolyPhen)
      ↓
Gene-wise Variant Extraction
(STAT6 and FCER1A)
      ↓
Variant Statistics
      ↓
Visualization Plots
      ↓
MultiQC HTML Report
      ↓
Final Analysis Report
```

---

## Pipeline Components

| Step                  | Tool Used        |
| --------------------- | ---------------- |
| SRA Download          | SRA Toolkit      |
| Quality Control       | FastQC           |
| Read Trimming         | fastp            |
| Alignment             | BWA-MEM          |
| BAM Processing        | SAMtools         |
| Variant Calling       | BCFtools         |
| Variant Annotation    | Ensembl VEP      |
| Functional Prediction | SIFT, PolyPhen   |
| Visualization         | GNUPlot / Python |
| HTML Report           | MultiQC          |

---

## Features

* Automated end-to-end NGS workflow
* Detection of variants in food allergy genes
* Gene-wise variant extraction
* SNP statistics generation
* Variant quality summary tables
* Functional annotation using VEP
* SIFT and PolyPhen prediction support
* Automatic visualization generation
* MultiQC HTML report generation
* Multi-sample analysis support
* Automated final report generation

---

## Input

* Paired-end FASTQ files downloaded from NCBI SRA.

Example SRA accession:

```text
SRR390728
SRR562646
SRR8534325
```

---

## Output Files

| File                     | Description            |
| ------------------------ | ---------------------- |
| variants.vcf             | Raw variant calls      |
| STAT6_variant_table.tsv  | STAT6 variants         |
| FCER1A_variant_table.tsv | FCER1A variants        |
| annotated_variants.txt   | Annotated variants     |
| snp_statistics.txt       | SNP distribution       |
| variant_counts.png       | Variant count plot     |
| multiqc_report.html      | Quality summary report |
| final_report.txt         | Final pipeline report  |

---

## Example Command

```bash
chmod +x pipeline.sh

./pipeline.sh
```

---

## Software Requirements

* SRA Toolkit
* FastQC
* fastp
* BWA
* SAMtools
* BCFtools
* Ensembl VEP
* MultiQC
* GNUPlot

---

## Future Improvements

* ANNOVAR-based annotation
* Machine-learning based variant prioritization
* Snakemake workflow implementation
* Docker containerization
* Comparative multi-sample analysis
* Interactive visualization dashboard

---

## Author

**Arushi Singh**

M.Tech Biotechnology

Thapar Institute of Engineering and Technology

--assembly GRCh38
