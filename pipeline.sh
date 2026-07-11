#!/bin/bash
#
# ============================================================================
#  Food Allergy Variant Analysis Pipeline v2.0
#  Target genes : STAT6 (chr12, GRCh38: 57,090,404-57,137,139)
#                  FCER1A (chr1, GRCh38: 159,278,575-159,313,224)
#  (5 kb flank added around each gene body to capture UTR / promoter
#   variants - verify against Ensembl/UCSC before publication use)
#
#  Reference    : Full GRCh38 primary assembly (no more custom mini-FASTA)
#  Purpose      : WGS/WES short-read alignment, duplicate marking,
#                 MAPQ filtering, variant calling, hard-filtering,
#                 VEP annotation, gene-region extraction, QC + reporting
# ============================================================================
#
# USAGE:
#   ./food_allergy_pipeline_v2.sh
#
# REQUIRED INPUT FILES (same directory as this script, or set paths below):
#   samples.txt                 - one SRA accession (or sample ID) per line
#   GRCh38_full_analysis_set.fa - full GRCh38 reference FASTA (see
#                                  prepare_reference() for the download URL)
#   target_genes.bed            - BED file with STAT6 / FCER1A coordinates
#                                  (auto-generated if missing - see below)
#
# DEPENDENCIES (must be on PATH):
#   bwa, samtools (>=1.15, for `samtools markdup`), bcftools (>=1.15),
#   fastqc, fastp, prefetch, fasterq-dump, vdb-dump (sra-tools),
#   gnuplot, multiqc, python3 (for plotting stats), tabix / bgzip (htslib)
#   vep (optional - annotation step is skipped with a warning if absent)
#   GNU parallel (optional - falls back to serial processing if absent)
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# 0. CONFIGURATION
# ---------------------------------------------------------------------------
# NOTE: keeping every tunable in one block makes the pipeline easy to adapt
# to a new project without hunting through the code (a common source of
# copy-paste bugs in the original script, e.g. the hard-coded -X 100000).

REF_DIR="reference"
REF="${REF_DIR}/GRCh38_full_analysis_set.fa"          # full genome, NOT a gene-slice
REF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
GENE_BED="target_genes.bed"                            # STAT6 + FCER1A regions
SAMPLES_FILE="samples.txt"

THREADS="$(nproc 2>/dev/null || echo 4)"               # auto-detect CPU count
MAX_PARALLEL_SAMPLES=2                                  # samples processed concurrently
                                                         # (keep low: bwa/mpileup are
                                                         #  already multi-threaded per sample)

# Variant filtering thresholds (GATK/bcftools community consensus for
# germline short-variant hard-filtering when a joint-genotyping/VQSR
# workflow is not used - see references section of the report)
MIN_QUAL=30
MIN_DP=10
MIN_MQ=30
MIN_MAPQ_READS=30                                       # read-level MAPQ filter (samtools view -q)

RESULTS_DIR="results"
LOG_DIR="logs"
CHECKPOINT_DIR=".checkpoints"                            # marks completed steps -> resumability
LOG_FILE="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

VEP_CACHE_DIR="${HOME}/.vep"                             # offline VEP cache location
VEP_ASSEMBLY="GRCh38"
VEP_SPECIES="homo_sapiens"

mkdir -p "$REF_DIR" "$RESULTS_DIR" "$LOG_DIR" "$CHECKPOINT_DIR"

# ---------------------------------------------------------------------------
# Logging + error handling
# ---------------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" | tee -a "$LOG_FILE" >&2
}

# Global trap only guards the main/setup steps. Per-sample failures are
# isolated inside process_sample() with their own error handling (see
# below) so that one bad sample does not kill the whole cohort - this was
# a real risk in the original script because `set -e` + a single global
# trap aborts the entire run on the first error in any sample.
on_error() {
    local exit_code=$?
    log "ERROR: Pipeline failed at line $1 (exit code $exit_code)"
    exit "$exit_code"
}
trap 'on_error $LINENO' ERR

# Checkpoint helpers: idempotent "has this step already run for this
# sample" markers, so the whole pipeline is safely re-runnable/resumable.
ckpt_path() { echo "${CHECKPOINT_DIR}/${1}.done"; }
is_done()   { [ -f "$(ckpt_path "$1")" ]; }
mark_done() { touch "$(ckpt_path "$1")"; }

# ---------------------------------------------------------------------------
# 1. Dependency checks
# ---------------------------------------------------------------------------
check_dependencies() {
    log "Checking required tools..."
    local required_tools=(bwa samtools bcftools fastqc fastp prefetch fasterq-dump vdb-dump gnuplot multiqc bgzip tabix python3)
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
        warn "VEP not found - annotation step will be SKIPPED. Install Ensembl VEP to enable SIFT/PolyPhen/CADD/ClinVar/gnomAD integration."
    fi

    if command -v parallel >/dev/null 2>&1; then
        PARALLEL_AVAILABLE=1
    else
        PARALLEL_AVAILABLE=0
        warn "GNU parallel not found - samples will be processed serially."
    fi

    # samtools markdup requires samtools >= 1.10 (name-sort + fixmate pipeline)
    local st_version
    st_version=$(samtools --version | head -1 | awk '{print $2}')
    log "samtools version: $st_version"

    log "Dependency check complete. THREADS=${THREADS}"
}

# ---------------------------------------------------------------------------
# 2. Reference genome: download, index, and build the gene-region BED
# ---------------------------------------------------------------------------
prepare_reference() {
    if [ ! -f "$REF" ]; then
        warn "Reference '$REF' not found."
        log "Download the full GRCh38 reference (recommended: the 1000 Genomes"
        log "'GRCh38_full_analysis_set_plus_decoy_hla.fa' build used by most"
        log "modern germline pipelines) with, e.g.:"
        log "    wget -O '${REF}' '${REF_URL}'"
        log "A whole-genome FASTA is large (~3 GB) and must be downloaded once"
        log "outside of this script (network access / storage permitting)."
        exit 1
    fi

    if [ ! -f "${REF}.bwt" ]; then
        log "Indexing reference with bwa index (this can take 1-2 h for full GRCh38)..."
        bwa index "$REF"
    else
        log "BWA index already exists. Skipping."
    fi

    if [ ! -f "${REF}.fai" ]; then
        log "Indexing reference with samtools faidx..."
        samtools faidx "$REF"
    else
        log "samtools .fai index already exists. Skipping."
    fi

    if [ ! -f "${REF%.fa}.dict" ] && [ ! -f "${REF%.fasta}.dict" ]; then
        log "Creating sequence dictionary (samtools dict)..."
        samtools dict "$REF" > "${REF%.*}.dict"
    fi

    # Gene-region BED (0-based, half-open, as required by bcftools/bedtools).
    # Coordinates below are GRCh38 (verified via Ensembl/GeneCards/NCBI,
    # July 2026) with a 5 kb flank to capture promoter/UTR/regulatory
    # variants. ALWAYS re-verify against Ensembl/UCSC before using for a
    # thesis or publication, since gene models are periodically revised.
    if [ ! -f "$GENE_BED" ]; then
        log "Creating target_genes.bed (STAT6, FCER1A; GRCh38 + 5kb flank)..."
        cat > "$GENE_BED" <<'BEDEOF'
chr12	57090404	57137139	STAT6
chr1	159278575	159313224	FCER1A
BEDEOF
        log "NOTE: if your reference FASTA uses RefSeq-style contig names"
        log "(e.g. NC_000012.12) instead of 'chr12', edit target_genes.bed"
        log "accordingly, or use 'bcftools annotate --rename-chrs' beforehand."
    fi
}

# ---------------------------------------------------------------------------
# 3. Per-sample: download + format detection
# ---------------------------------------------------------------------------
download_sample() {
    local SAMPLE="$1"
    local step="download_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] Download already complete (checkpoint)."; return 0; fi

    log "[$SAMPLE] Detecting SRA data type..."
    local FORMAT
    FORMAT=$(vdb-dump "$SAMPLE" --info 2>/dev/null | awk -F': ' '/FMT/{print $2}') || FORMAT="UNKNOWN"
    log "[$SAMPLE] Detected format: ${FORMAT:-UNKNOWN}"

    case "$FORMAT" in
        FASTQ|SRA)
            if [ ! -f "${SAMPLE}_1.fastq.gz" ] && [ ! -f "${SAMPLE}_1.fastq" ]; then
                prefetch "$SAMPLE"
                # NOTE: the original script capped downloads at -X 100000
                # reads. That silently analyzes only a subset of the run
                # and will under-report both coverage and rare variants.
                # Full-depth download by default; --gzip saves disk space.
                fasterq-dump "$SAMPLE" -e "$THREADS" -p --split-files
                gzip -f "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq"
            fi
            SAMPLE_TYPE="FASTQ"
            ;;
        BAM)
            if [ ! -f "${SAMPLE}.bam" ]; then
                prefetch "$SAMPLE"
                sam-dump "$SAMPLE" | samtools view -bS - > "${SAMPLE}.bam"
            fi
            SAMPLE_TYPE="BAM"
            ;;
        *)
            warn "[$SAMPLE] Unsupported/undetected format ('$FORMAT'). Skipping sample."
            return 1
            ;;
    esac
    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 4. QC + trimming (FastQC pre -> fastp -> FastQC post)
# ---------------------------------------------------------------------------
qc_and_trim() {
    local SAMPLE="$1"
    local step="qc_trim_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] QC/trim already complete (checkpoint)."; return 0; fi
    [ "$SAMPLE_TYPE" = "BAM" ] && { log "[$SAMPLE] BAM input - skipping FASTQ QC/trim."; mark_done "$step"; return 0; }

    local R1="${SAMPLE}_1.fastq.gz" R2="${SAMPLE}_2.fastq.gz"

    log "[$SAMPLE] FastQC on RAW reads..."
    mkdir -p "${RESULTS_DIR}/${SAMPLE}/fastqc_raw"
    fastqc -t "$THREADS" -o "${RESULTS_DIR}/${SAMPLE}/fastqc_raw" "$R1" "$R2"

    log "[$SAMPLE] Trimming with fastp..."
    fastp \
        -i "$R1" -I "$R2" \
        -o "${SAMPLE}_1_clean.fastq.gz" -O "${SAMPLE}_2_clean.fastq.gz" \
        --cut_right --cut_right_mean_quality 20 \
        --length_required 30 \
        --detect_adapter_for_pe \
        -w "$THREADS" \
        -h "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_fastp.html" \
        -j "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_fastp.json"

    log "[$SAMPLE] FastQC on TRIMMED reads (post-trim check requested by the user)..."
    mkdir -p "${RESULTS_DIR}/${SAMPLE}/fastqc_trimmed"
    fastqc -t "$THREADS" -o "${RESULTS_DIR}/${SAMPLE}/fastqc_trimmed" \
        "${SAMPLE}_1_clean.fastq.gz" "${SAMPLE}_2_clean.fastq.gz"

    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 5. Alignment, sort, duplicate marking, MAPQ filter, index
# ---------------------------------------------------------------------------
align_and_process_bam() {
    local SAMPLE="$1"
    local step="align_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] Alignment/BAM processing already complete."; return 0; fi

    if [ "$SAMPLE_TYPE" = "FASTQ" ]; then
        log "[$SAMPLE] Aligning with BWA-MEM..."
        # -R sets a proper read group; bcftools/VEP/Picard-style tools
        # generally expect one for multi-sample cohorts.
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "$REF" "${SAMPLE}_1_clean.fastq.gz" "${SAMPLE}_2_clean.fastq.gz" \
            | samtools sort -@ "$THREADS" -o "${SAMPLE}_sorted.bam" -
    else
        log "[$SAMPLE] BAM input already aligned - sorting by coordinate..."
        samtools sort -@ "$THREADS" -o "${SAMPLE}_sorted.bam" "${SAMPLE}.bam"
    fi

    # --- Duplicate marking -------------------------------------------------
    # The original script never removed PCR duplicates. Unmarked
    # duplicates inflate depth (DP) and can produce false-positive
    # variant calls, especially at PCR-amplified allergy-panel loci.
    # `samtools markdup` requires name-sorted + fixmate input first.
    log "[$SAMPLE] Marking PCR/optical duplicates (samtools markdup)..."
    samtools sort -@ "$THREADS" -n -o "${SAMPLE}_namesorted.bam" "${SAMPLE}_sorted.bam"
    samtools fixmate -@ "$THREADS" -m "${SAMPLE}_namesorted.bam" "${SAMPLE}_fixmate.bam"
    samtools sort -@ "$THREADS" -o "${SAMPLE}_fixmate_sorted.bam" "${SAMPLE}_fixmate.bam"
    samtools markdup -@ "$THREADS" -s \
        "${SAMPLE}_fixmate_sorted.bam" "${SAMPLE}_dedup.bam" \
        2> "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_markdup_stats.txt"
    rm -f "${SAMPLE}_namesorted.bam" "${SAMPLE}_fixmate.bam" "${SAMPLE}_fixmate_sorted.bam"
    # (Alternative: `picard MarkDuplicates I=... O=... M=metrics.txt`
    #  produces a richer per-library metrics report if Picard is available.)

    # --- MAPQ filtering ------------------------------------------------------
    # Removes multi-mapping / ambiguously-placed reads (MAPQ < 30) before
    # variant calling, reducing false calls from mis-mapped reads in
    # paralogous or repetitive regions.
    log "[$SAMPLE] Filtering reads with MAPQ >= ${MIN_MAPQ_READS}..."
    samtools view -@ "$THREADS" -b -q "$MIN_MAPQ_READS" -F 0x400 \
        "${SAMPLE}_dedup.bam" > "${SAMPLE}_final.bam"
        # -F 0x400 additionally drops reads already flagged as duplicates

    samtools index -@ "$THREADS" "${SAMPLE}_final.bam"

    # QC metrics on the final analysis-ready BAM
    samtools flagstat "${SAMPLE}_final.bam" > "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_flagstat.txt"
    samtools idxstats "${SAMPLE}_final.bam" > "${RESULTS_DIR}/${SAMPLE}/${SAMPLE}_idxstats.txt"

    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 6. Variant calling + hard filtering
# ---------------------------------------------------------------------------
call_variants() {
    local SAMPLE="$1"
    local step="call_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] Variant calling already complete."; return 0; fi

    log "[$SAMPLE] Calling variants (bcftools mpileup | bcftools call)..."
    # -a AD,DP,MQ : retain per-allele depth, total depth, mapping quality
    #               so they are available for hard-filtering downstream.
    # -Ou / -Oz   : uncompressed BCF between piped steps for speed, then
    #               compressed+indexed VCF for storage/downstream tools.
    bcftools mpileup \
        -f "$REF" \
        -a FORMAT/AD,FORMAT/DP,INFO/MQ \
        --max-depth 250 \
        --threads "$THREADS" \
        -Ou "${SAMPLE}_final.bam" | \
    bcftools call \
        -mv \
        --threads "$THREADS" \
        -Oz -o "${SAMPLE}_raw_variants.vcf.gz"
    tabix -p vcf "${SAMPLE}_raw_variants.vcf.gz"

    # --- Hard filtering ------------------------------------------------------
    # The original pipeline never filtered variants at all, which mixes
    # low-confidence/noise calls with real allergy-associated variants.
    # QUAL>30, DP>10, MQ>30 are conventional community thresholds for a
    # single-sample bcftools workflow without VQSR/joint genotyping.
    log "[$SAMPLE] Applying hard filters (QUAL>${MIN_QUAL}, DP>${MIN_DP}, MQ>${MIN_MQ})..."
    bcftools filter \
        -e "QUAL<${MIN_QUAL} || INFO/DP<${MIN_DP} || INFO/MQ<${MIN_MQ}" \
        -s LOWQUAL -Oz -o "${SAMPLE}_filtered.vcf.gz" \
        "${SAMPLE}_raw_variants.vcf.gz"
    tabix -p vcf "${SAMPLE}_filtered.vcf.gz"

    # Keep only variants that PASS for the downstream annotation/report,
    # but retain the flagged (LOWQUAL) file too for transparency/QC.
    bcftools view -f PASS -Oz -o "${SAMPLE}_pass.vcf.gz" "${SAMPLE}_filtered.vcf.gz"
    tabix -p vcf "${SAMPLE}_pass.vcf.gz"

    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 7. VEP annotation (SIFT / PolyPhen-2 / CADD / ClinVar / gnomAD)
# ---------------------------------------------------------------------------
annotate_variants() {
    local SAMPLE="$1"
    local step="annotate_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] Annotation already complete."; return 0; fi
    if [ "${VEP_AVAILABLE:-0}" -ne 1 ]; then
        warn "[$SAMPLE] Skipping VEP annotation (VEP not installed)."
        cp "${SAMPLE}_pass.vcf.gz" "${SAMPLE}_annotated.vcf.gz"
        mark_done "$step"
        return 0
    fi

    log "[$SAMPLE] Running Ensembl VEP..."
    # --sift b / --polyphen b : return both prediction + score
    # --plugin CADD            : requires locally downloaded CADD score files
    # --custom (ClinVar)       : requires a local bgzipped/tabixed ClinVar VCF
    # --custom (gnomAD)        : requires a local bgzipped/tabixed gnomAD VCF
    # See the report's "Explanation" section for how to obtain and wire in
    # each of these resource files; they are NOT bundled with VEP itself.
    vep \
        --input_file "${SAMPLE}_pass.vcf.gz" \
        --output_file "${SAMPLE}_annotated.vcf" \
        --vcf --force_overwrite --offline --cache \
        --dir_cache "$VEP_CACHE_DIR" \
        --assembly "$VEP_ASSEMBLY" \
        --species "$VEP_SPECIES" \
        --fasta "$REF" \
        --sift b --polyphen b \
        --symbol --canonical --numbers \
        $( [ -d "${VEP_CACHE_DIR}/Plugins" ] && echo "--plugin CADD,${VEP_CACHE_DIR}/CADD/whole_genome_SNVs.tsv.gz" ) \
        $( [ -f "${VEP_CACHE_DIR}/clinvar/clinvar.vcf.gz" ] && echo "--custom ${VEP_CACHE_DIR}/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT" ) \
        $( [ -f "${VEP_CACHE_DIR}/gnomad/gnomad.genomes.vcf.gz" ] && echo "--custom ${VEP_CACHE_DIR}/gnomad/gnomad.genomes.vcf.gz,gnomAD,vcf,exact,0,AF" ) \
        2>> "$LOG_FILE" || warn "[$SAMPLE] VEP annotation encountered errors (continuing with unannotated PASS VCF)."

    if [ -f "${SAMPLE}_annotated.vcf" ]; then
        bgzip -f "${SAMPLE}_annotated.vcf"
        tabix -p vcf "${SAMPLE}_annotated.vcf.gz"
    else
        cp "${SAMPLE}_pass.vcf.gz" "${SAMPLE}_annotated.vcf.gz"
    fi
    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 8. Gene-region extraction, statistics, plots, per-sample report
# ---------------------------------------------------------------------------
summarize_sample() {
    local SAMPLE="$1"
    local step="summarize_${SAMPLE}"
    if is_done "$step"; then log "[$SAMPLE] Summary already complete."; return 0; fi

    local OUT="${RESULTS_DIR}/${SAMPLE}"
    mkdir -p "$OUT"

    # --- Region extraction ---------------------------------------------------
    # The original script grep'd on a *whole-chromosome* accession token
    # (e.g. NC_000012), which pulls in every chromosome-12 variant, not
    # just the STAT6 gene. Using `bcftools view -R <BED>` on an indexed
    # VCF restricts to the exact gene +/-flank coordinates and is robust
    # to header/formatting differences that break plain-text grep.
    log "[$SAMPLE] Extracting STAT6 / FCER1A variants via bcftools regions..."
    bcftools view -R "$GENE_BED" -Oz -o "${OUT}/${SAMPLE}_target_genes.vcf.gz" \
        "${SAMPLE}_annotated.vcf.gz"
    tabix -p vcf "${OUT}/${SAMPLE}_target_genes.vcf.gz"

    awk -F'\t' '$4=="STAT6"' "$GENE_BED" > /tmp/stat6.bed
    awk -F'\t' '$4=="FCER1A"' "$GENE_BED" > /tmp/fcer1a.bed
    bcftools view -R /tmp/stat6.bed -Ov "${OUT}/${SAMPLE}_target_genes.vcf.gz" \
        > "${OUT}/${SAMPLE}_STAT6.vcf"
    bcftools view -R /tmp/fcer1a.bed -Ov "${OUT}/${SAMPLE}_target_genes.vcf.gz" \
        > "${OUT}/${SAMPLE}_FCER1A.vcf"

    # Tab-separated summary tables (beginner-friendly, spreadsheet-ready)
    {
        echo -e "Chromosome\tPosition\tREF\tALT\tQUAL\tDP\tMQ"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/MQ\n' \
            "${OUT}/${SAMPLE}_STAT6.vcf" 2>/dev/null || true
    } > "${OUT}/${SAMPLE}_STAT6_table.tsv"

    {
        echo -e "Chromosome\tPosition\tREF\tALT\tQUAL\tDP\tMQ"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/MQ\n' \
            "${OUT}/${SAMPLE}_FCER1A.vcf" 2>/dev/null || true
    } > "${OUT}/${SAMPLE}_FCER1A_table.tsv"

    # --- bcftools stats: SNP/INDEL counts + Ts/Tv ratio ----------------------
    # `bcftools stats` computes Ts/Tv, SNP/indel counts, quality
    # distributions etc. in one pass - far more reliable than hand-rolled
    # awk counting, and it's what MultiQC parses automatically.
    log "[$SAMPLE] Computing variant statistics (bcftools stats)..."
    bcftools stats "${OUT}/${SAMPLE}_target_genes.vcf.gz" > "${OUT}/${SAMPLE}_bcftools_stats.txt"

    local STAT6_COUNT FCER1A_COUNT TOTAL TSTV
    STAT6_COUNT=$(bcftools view -H "${OUT}/${SAMPLE}_STAT6.vcf" | wc -l)
    FCER1A_COUNT=$(bcftools view -H "${OUT}/${SAMPLE}_FCER1A.vcf" | wc -l)
    TOTAL=$(bcftools view -H "${OUT}/${SAMPLE}_target_genes.vcf.gz" | wc -l)
    TSTV=$(grep "^TSTV" "${OUT}/${SAMPLE}_bcftools_stats.txt" | awk -F'\t' '{print $5}' | head -1)

    # --- Publication-quality plots (Python/matplotlib instead of raw
    #     gnuplot) - variant counts per gene, SNP vs INDEL, quality
    #     distribution, and variant density along each gene. -----------------
    python3 "$(dirname "$0")/plot_variant_summary.py" \
        --sample "$SAMPLE" \
        --stats "${OUT}/${SAMPLE}_bcftools_stats.txt" \
        --stat6-vcf "${OUT}/${SAMPLE}_STAT6.vcf" \
        --fcer1a-vcf "${OUT}/${SAMPLE}_FCER1A.vcf" \
        --outdir "$OUT" \
        || warn "[$SAMPLE] Plot generation failed (continuing)."

    # --- Per-sample final report ---------------------------------------------
    {
        echo "Food Allergy Variant Analysis Report"
        echo "====================================="
        echo ""
        echo "Sample ID: $SAMPLE"
        echo "Reference: GRCh38 (full primary assembly)"
        echo "Total PASS variants in target regions: $TOTAL"
        echo "STAT6 variants: $STAT6_COUNT"
        echo "FCER1A variants: $FCER1A_COUNT"
        echo "Transition/Transversion (Ts/Tv) ratio: ${TSTV:-NA}"
        echo ""
        echo "Filters applied: QUAL>${MIN_QUAL}, DP>${MIN_DP}, MQ>${MIN_MQ}, read MAPQ>=${MIN_MAPQ_READS}, duplicates removed"
        echo ""
        echo "See ${SAMPLE}_bcftools_stats.txt for full SNP/INDEL/Ts-Tv breakdown."
    } > "${OUT}/${SAMPLE}_final_report.txt"

    log "[$SAMPLE] ===== SUMMARY ====="
    cat "${OUT}/${SAMPLE}_final_report.txt" | tee -a "$LOG_FILE"

    mark_done "$step"
}

# ---------------------------------------------------------------------------
# 9. Cleanup (keep checkpointed outputs, remove large intermediates only)
# ---------------------------------------------------------------------------
cleanup_sample() {
    local SAMPLE="$1"
    log "[$SAMPLE] Cleaning up large intermediates (raw/trimmed FASTQ, SRA cache)..."
    rm -f "${SAMPLE}_1.fastq.gz" "${SAMPLE}_2.fastq.gz" \
          "${SAMPLE}_1_clean.fastq.gz" "${SAMPLE}_2_clean.fastq.gz" \
          "${SAMPLE}_dedup.bam" "${SAMPLE}.bam"
    rm -rf "${SAMPLE}"   # sra-tools prefetch cache directory
    # NOTE: ${SAMPLE}_final.bam / _final.bam.bai and all VCFs are KEPT -
    # they are the checkpointed, reusable analysis-ready outputs.
}

# ---------------------------------------------------------------------------
# 10. Orchestration for one sample (isolated error handling)
# ---------------------------------------------------------------------------
process_sample() {
    local SAMPLE="$1"
    log "========================================="
    log "Processing Sample: $SAMPLE"
    log "========================================="

    # Run this sample's steps in a subshell with its own error trap so a
    # failure here is logged and skipped, NOT fatal to the whole cohort.
    (
        set -e
        SAMPLE_TYPE="FASTQ"
        download_sample "$SAMPLE"
        qc_and_trim "$SAMPLE"
        align_and_process_bam "$SAMPLE"
        call_variants "$SAMPLE"
        annotate_variants "$SAMPLE"
        summarize_sample "$SAMPLE"
        cleanup_sample "$SAMPLE"
    ) || { warn "[$SAMPLE] FAILED - see $LOG_FILE for details. Continuing with next sample."; return 0; }

    log "Finished processing sample $SAMPLE."
}
export -f process_sample download_sample qc_and_trim align_and_process_bam \
           call_variants annotate_variants summarize_sample cleanup_sample \
           log warn is_done mark_done ckpt_path
export LOG_FILE REF GENE_BED THREADS MIN_QUAL MIN_DP MIN_MQ MIN_MAPQ_READS \
       RESULTS_DIR CHECKPOINT_DIR VEP_AVAILABLE VEP_CACHE_DIR VEP_ASSEMBLY VEP_SPECIES

# ---------------------------------------------------------------------------
# 11. Main
# ---------------------------------------------------------------------------
main() {
    log "========================================="
    log " FOOD ALLERGY NGS PIPELINE v2.0 STARTED"
    log "========================================="

    check_dependencies
    prepare_reference

    if [ ! -f "$SAMPLES_FILE" ]; then log "Samples file '$SAMPLES_FILE' not found!"; exit 1; fi
    if [ ! -s "$SAMPLES_FILE" ]; then log "Samples file '$SAMPLES_FILE' is empty!"; exit 1; fi

    if [ "${PARALLEL_AVAILABLE:-0}" -eq 1 ]; then
        log "Processing samples in parallel (max ${MAX_PARALLEL_SAMPLES} concurrent)..."
        parallel -j "$MAX_PARALLEL_SAMPLES" process_sample :::: "$SAMPLES_FILE"
    else
        while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]; do
            [ -z "$SAMPLE" ] && continue
            process_sample "$SAMPLE"
        done < "$SAMPLES_FILE"
    fi

    log "Running MultiQC aggregate report across all samples..."
    multiqc "$RESULTS_DIR" -o "${RESULTS_DIR}/multiqc_output"

    log "========================================="
    log " PIPELINE COMPLETED"
    log " Per-sample outputs: ${RESULTS_DIR}/<SAMPLE>/"
    log " Cohort QC summary : ${RESULTS_DIR}/multiqc_output/"
    log "========================================="
}

main "$@"
