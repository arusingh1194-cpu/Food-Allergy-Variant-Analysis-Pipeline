#!/bin/bash
#
# ============================================================
# Food Allergy Variant Analysis Pipeline
# Target genes: STAT6 (chr12) and FCER1A (chr1)
# Optimized for WSL / local laptop execution
# ============================================================
#
# Usage:
#   ./food_allergy_pipeline.sh
#
# Requirements:
#   samples.txt          - one SRA accession per line
#   STAT6_FCER1A.fasta   - reference FASTA covering both genes
#
# Dependencies (must be on PATH):
#   bwa, samtools, bcftools, fastqc, fastp,
#   prefetch, fasterq-dump (sra-tools), gnuplot, multiqc
#   vep (optional - annotation step is skipped if absent)
#
set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
REF="STAT6_FCER1A.fasta"
SAMPLES_FILE="samples.txt"
THREADS=4
LOG_FILE="pipeline_$(date +%Y%m%d_%H%M%S).log"

# Chromosome accession tokens used to split variants by gene.
# NOTE: adjust these if your reference uses different contig
# naming (e.g. "chr1" / "chr12" instead of RefSeq accessions).
STAT6_CHROM_TOKEN="NC_000012"
FCER1A_CHROM_TOKEN="NC_000001"

# ------------------------------------------------------------
# Logging helper
# ------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

on_error() {
    local exit_code=$?
    log "ERROR: Pipeline failed at line $1 with exit code $exit_code"
    exit "$exit_code"
}
trap 'on_error $LINENO' ERR

# ------------------------------------------------------------
# Step 0: Dependency checks
# ------------------------------------------------------------
check_dependencies() {
    log "Checking required tools..."
    local required_tools=(bwa samtools bcftools fastqc fastp prefetch fasterq-dump gnuplot multiqc)
    local missing=0

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "MISSING required tool: $tool"
            missing=1
        fi
    done

    if [ "$missing" -eq 1 ]; then
        log "One or more required tools are missing. Aborting."
        exit 1
    fi

    if command -v vep >/dev/null 2>&1; then
        VEP_AVAILABLE=1
        log "VEP found - annotation step will run."
    else
        VEP_AVAILABLE=0
        log "VEP not found - annotation step will be skipped."
    fi

    log "All required tools are available."
}

# ------------------------------------------------------------
# Step 1: Reference genome checks and indexing
# ------------------------------------------------------------
prepare_reference() {
    if [ ! -f "$REF" ]; then
        log "Reference genome '$REF' not found!"
        exit 1
    fi

    if [ ! -f "${REF}.bwt" ]; then
        log "Indexing reference genome with bwa..."
        bwa index "$REF"
    else
        log "BWA index already exists. Skipping."
    fi

    if [ ! -f "${REF}.fai" ]; then
        log "Indexing reference genome with samtools faidx..."
        samtools faidx "$REF"
    else
        log "samtools .fai index already exists. Skipping."
    fi
}

# ------------------------------------------------------------
# Step 2: Per-sample processing
# ------------------------------------------------------------
process_sample() {
    local SAMPLE="$1"

    log "========================================="
    log "Processing Sample: $SAMPLE"
    log "========================================="
########################################################
# DETECT DATA TYPE
########################################################

echo "Detecting input format..."

FORMAT=$(vdb-dump "$SAMPLE" --info 2>/dev/null | awk -F': ' '/FMT/{print $2}')

echo "Detected format: $FORMAT"

########################################################
# DOWNLOAD DATA
########################################################

if [[ "$FORMAT" == "FASTQ" || "$FORMAT" == "SRA" ]]; then

    echo "Raw sequencing data detected."

    if [[ ! -f "${SAMPLE}_1.fastq" && ! -f "${SAMPLE}_1.fastq.gz" ]]; then

        prefetch "$SAMPLE"

        fasterq-dump "$SAMPLE" \
            -e "$THREADS" \
            -p

    fi

    READ1="${SAMPLE}_1.fastq"
    READ2="${SAMPLE}_2.fastq"

elif [[ "$FORMAT" == "BAM" ]]; then

    echo "Aligned BAM dataset detected."

    if [[ ! -f "${SAMPLE}.bam" ]]; then

        prefetch "$SAMPLE"

        sam-dump "$SAMPLE" | samtools view -bS - > "${SAMPLE}.bam"

    fi

    cp "${SAMPLE}.bam" "${SAMPLE}_sorted.bam"

    samtools index "${SAMPLE}_sorted.bam"

    BAM_ALREADY_AVAILABLE=1

else

    echo "Unsupported format: $FORMAT"

    continue

fi
    # --- Download ---
    if [ ! -f "${SAMPLE}_1.fastq" ] || [ ! -f "${SAMPLE}_2.fastq" ]; then
        log "Downloading sample $SAMPLE..."
        prefetch "$SAMPLE"
        fasterq-dump "$SAMPLE" \
            -X 100000 \
            -e "$THREADS" \
            -p

        if [ ! -f "${SAMPLE}_1.fastq" ] || [ ! -f "${SAMPLE}_2.fastq" ]; then
            log "FASTQ download failed for $SAMPLE. Skipping sample."
            return 0
        fi
    else
        log "FASTQ files already present for $SAMPLE. Skipping download."
    fi
if [[ -z "${BAM_ALREADY_AVAILABLE:-}" ]]; then
    fastqc "$READ1" "$READ2"

    fastp \
        -i "$READ1" \
        -I "$READ2" \
        -o "${SAMPLE}_1_clean.fastq" \
        -O "${SAMPLE}_2_clean.fastq" \
        -w "$THREADS"

fi
    # --- FastQC (raw) ---
    if [ ! -f "${SAMPLE}_1_fastqc.html" ]; then
        log "Running FastQC on raw reads..."
        fastqc "${SAMPLE}_1.fastq"
        fastqc "${SAMPLE}_2.fastq"
    else
        log "FastQC results already exist. Skipping."
    fi

    # --- Trimming ---
    if [ ! -f "${SAMPLE}_1_clean.fastq" ]; then
        log "Trimming reads with fastp..."
        fastp \
            -i "${SAMPLE}_1.fastq" \
            -I "${SAMPLE}_2.fastq" \
            -o "${SAMPLE}_1_clean.fastq" \
            -O "${SAMPLE}_2_clean.fastq" \
            --cut_right \
            --cut_right_mean_quality 20 \
            --length_required 30 \
            -w "$THREADS" \
            -h "${SAMPLE}_fastp.html" \
            -j "${SAMPLE}_fastp.json"
    else
        log "Trimmed reads already exist. Skipping fastp."
    fi
########################################################
# ALIGNMENT
########################################################

if [[ -z "${BAM_ALREADY_AVAILABLE:-}" ]]; then

    bwa mem \
        -t "$THREADS" \
        "$REF" \
        "$READ1" \
        "$READ2" | \
    samtools sort \
        -@ "$THREADS" \
        -o "${SAMPLE}_sorted.bam"

    samtools index "${SAMPLE}_sorted.bam"

else

    echo "Alignment skipped (input already contains aligned BAM)."

fi

    # --- BAM index ---
    if [ ! -f "${SAMPLE}_sorted.bam.bai" ]; then
        log "Indexing BAM file..."
        samtools index "${SAMPLE}_sorted.bam"
    else
        log "BAM index already exists. Skipping."
    fi

    # --- Variant calling ---
    if [ ! -f "${SAMPLE}_variants.vcf" ]; then
        log "Calling variants with bcftools..."
        bcftools mpileup \
            -f "$REF" \
            "${SAMPLE}_sorted.bam" | \
        bcftools call \
            -mv \
            -Ov \
            -o "${SAMPLE}_variants.vcf"
    else
        log "Variant file already exists. Skipping variant calling."
    fi

    # --- Optional VEP annotation ---
    if [ "${VEP_AVAILABLE:-0}" -eq 1 ]; then
        if [ ! -f "${SAMPLE}_variants.vep.vcf" ]; then
            log "Running VEP annotation..."
            vep \
                -i "${SAMPLE}_variants.vcf" \
                -o "${SAMPLE}_variants.vep.vcf" \
                --fasta "$REF" \
                --vcf \
                --force_overwrite || log "VEP annotation failed for $SAMPLE (continuing)."
        else
            log "VEP annotation already exists. Skipping."
        fi
    fi

    # --- Gene-specific extraction ---
    log "Extracting gene-specific variants..."
    grep "$STAT6_CHROM_TOKEN" "${SAMPLE}_variants.vcf" > "${SAMPLE}_STAT6.txt" || true
    grep "$FCER1A_CHROM_TOKEN" "${SAMPLE}_variants.vcf" > "${SAMPLE}_FCER1A.txt" || true

    local STAT6_COUNT FCER1A_COUNT TOTAL
    STAT6_COUNT=$(grep -vc "^#" "${SAMPLE}_STAT6.txt" || true)
    FCER1A_COUNT=$(grep -vc "^#" "${SAMPLE}_FCER1A.txt" || true)
    TOTAL=$(grep -vc "^#" "${SAMPLE}_variants.vcf" || true)

    # --- Variant tables ---
    log "Building variant tables..."
    {
        echo -e "Chromosome\tPosition\tREF\tALT\tQUAL"
        grep -v "^#" "${SAMPLE}_STAT6.txt" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}'
    } > "${SAMPLE}_STAT6_table.tsv"

    {
        echo -e "Chromosome\tPosition\tREF\tALT\tQUAL"
        grep -v "^#" "${SAMPLE}_FCER1A.txt" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}'
    } > "${SAMPLE}_FCER1A_table.tsv"

    # --- SNP statistics ---
    log "Computing SNP statistics..."
    grep -v "^#" "${SAMPLE}_variants.vcf" | \
        awk '{print $4">"$5}' | \
        sort | uniq -c > "${SAMPLE}_snp_statistics.txt"

    # --- Plot generation (per-sample counts file to avoid overwrite bug) ---
    log "Generating variant count plot..."
    local COUNTS_FILE="${SAMPLE}_counts.txt"
    {
        echo "STAT6 $STAT6_COUNT"
        echo "FCER1A $FCER1A_COUNT"
    } > "$COUNTS_FILE"

    gnuplot << EOF
set terminal png
set output "${SAMPLE}_variant_counts.png"
set style data histograms
set style fill solid
set xlabel "Genes"
set ylabel "Variant Count"
plot '${COUNTS_FILE}' using 2:xtic(1) title 'Variants'
EOF

    # --- Final report ---
    log "Writing final report..."
    {
        echo "Food Allergy Variant Analysis Report"
        echo "==================================="
        echo ""
        echo "Sample ID: $SAMPLE"
        echo "Total Variants: $TOTAL"
        echo "STAT6 Variants: $STAT6_COUNT"
        echo "FCER1A Variants: $FCER1A_COUNT"
        echo ""
        echo "SNP Distribution"
        cat "${SAMPLE}_snp_statistics.txt"
    } > "${SAMPLE}_final_report.txt"

    # --- Display results ---
    log "===== REPORT (${SAMPLE}) ====="
    cat "${SAMPLE}_final_report.txt" | tee -a "$LOG_FILE"

    log "===== STAT6 VARIANTS (head) ====="
    head "${SAMPLE}_STAT6_table.tsv" | tee -a "$LOG_FILE"

    log "===== FCER1A VARIANTS (head) ====="
    head "${SAMPLE}_FCER1A_table.tsv" | tee -a "$LOG_FILE"

    # --- Cleanup: remove only large/temporary intermediates ---
    # Raw FASTQs and SRA cache are safe to remove since downstream
    # trimmed/aligned/called files already exist and are checkpointed.
    log "Cleaning up temporary files for $SAMPLE..."
    rm -f "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq"
    rm -rf "${SAMPLE}"   # sra-tools prefetch cache directory

    log "Finished processing sample $SAMPLE."
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
main() {
    log "========================================="
    log " FOOD ALLERGY NGS PIPELINE STARTED"
    log "========================================="

    check_dependencies
    prepare_reference

    if [ ! -f "$SAMPLES_FILE" ]; then
        log "Samples file '$SAMPLES_FILE' not found!"
        exit 1
    fi

    if [ ! -s "$SAMPLES_FILE" ]; then
        log "Samples file '$SAMPLES_FILE' is empty!"
        exit 1
    fi

    while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]; do
        [ -z "$SAMPLE" ] && continue
        process_sample "$SAMPLE"
    done < "$SAMPLES_FILE"

    log "Running MultiQC aggregate report..."
    multiqc . -o multiqc_output

    log "========================================="
    log " PIPELINE COMPLETED SUCCESSFULLY"
    log "========================================="
}

main "$@"
