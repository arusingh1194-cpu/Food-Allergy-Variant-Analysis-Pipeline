#!/bin/bash

########################################################
# FOOD ALLERGY VARIANT ANALYSIS PIPELINE
# STAT6 and FCER1A
########################################################

echo "========================================="
echo " FOOD ALLERGY NGS PIPELINE STARTED"
echo "========================================="

SAMPLE=SRR562646
REF=STAT6_FCER1A.fasta
THREADS=4

########################################################
# STEP 1: DOWNLOAD DATA
########################################################

if [ ! -f "${SAMPLE}_1.fastq" ]; then

    echo "Downloading SRA sample..."

    prefetch $SAMPLE
    fasterq-dump $SAMPLE -e $THREADS

else
    echo "FASTQ files already exist."
fi

########################################################
# STEP 2: FASTQC
########################################################

if [ ! -f "${SAMPLE}_1_fastqc.html" ]; then

    fastqc ${SAMPLE}_1.fastq
    fastqc ${SAMPLE}_2.fastq

fi

########################################################
# STEP 3: TRIMMING
########################################################

if [ ! -f "${SAMPLE}_1_clean.fastq" ]; then

    fastp \
    -i ${SAMPLE}_1.fastq \
    -I ${SAMPLE}_2.fastq \
    -o ${SAMPLE}_1_clean.fastq \
    -O ${SAMPLE}_2_clean.fastq \
    --cut_right \
    --cut_right_mean_quality 20 \
    --length_required 30 \
    -h fastp_report.html \
    -j fastp_report.json

fi

########################################################
# STEP 4: FASTQC AFTER CLEANING
########################################################

if [ ! -f "${SAMPLE}_1_clean_fastqc.html" ]; then

    fastqc ${SAMPLE}_1_clean.fastq
    fastqc ${SAMPLE}_2_clean.fastq

fi

########################################################
# STEP 5: REFERENCE INDEX
########################################################

if [ ! -f "${REF}.bwt" ]; then

    bwa index $REF
    samtools faidx $REF

fi

########################################################
# STEP 6: ALIGNMENT + SORTING
########################################################

if [ ! -f "aligned_sorted.bam" ]; then

    bwa mem -t $THREADS \
    $REF \
    ${SAMPLE}_1_clean.fastq \
    ${SAMPLE}_2_clean.fastq | \
    samtools sort -@ $THREADS \
    -o aligned_sorted.bam

fi

########################################################
# STEP 7: BAM INDEX
########################################################

if [ ! -f "aligned_sorted.bam.bai" ]; then

    samtools index aligned_sorted.bam

fi

########################################################
# STEP 8: VARIANT CALLING
########################################################

if [ ! -f "variants.vcf" ]; then

    bcftools mpileup \
    -f $REF \
    aligned_sorted.bam | \
    bcftools call -mv \
    -Ov \
    -o variants.vcf

fi

########################################################
# STEP 9: TOTAL VARIANT COUNT
########################################################

TOTAL=$(grep -v "^#" variants.vcf | wc -l)

echo "Total variants = $TOTAL"

########################################################
# STEP 10: GENE-WISE EXTRACTION
########################################################

grep "NC_000012" variants.vcf > STAT6_variants.txt
grep "NC_000001" variants.vcf > FCER1A_variants.txt

STAT6_COUNT=$(grep -v "^#" STAT6_variants.txt | wc -l)
FCER1A_COUNT=$(grep -v "^#" FCER1A_variants.txt | wc -l)

########################################################
# STEP 11: VARIANT TABLES
########################################################

echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \
> STAT6_variant_table.tsv

grep -v "^#" STAT6_variants.txt | \
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
>> STAT6_variant_table.tsv


echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \
> FCER1A_variant_table.tsv

grep -v "^#" FCER1A_variants.txt | \
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
>> FCER1A_variant_table.tsv

########################################################
# STEP 12: SNP TYPE STATISTICS
########################################################

grep -v "^#" variants.vcf | \
awk '{print $4">"$5}' | \
sort | uniq -c > snp_statistics.txt

########################################################
# STEP 13: PIE CHART DATA
########################################################

echo "STAT6 $STAT6_COUNT" > piechart_data.txt
echo "FCER1A $FCER1A_COUNT" >> piechart_data.txt

########################################################
# STEP 14: FINAL REPORT
########################################################

echo "Food Allergy Variant Analysis Report" \
> final_report.txt

echo "===================================" \
>> final_report.txt

echo "" >> final_report.txt

echo "Sample ID: $SAMPLE" \
>> final_report.txt

echo "Total Variants: $TOTAL" \
>> final_report.txt

echo "STAT6 Variants: $STAT6_COUNT" \
>> final_report.txt

echo "FCER1A Variants: $FCER1A_COUNT" \
>> final_report.txt

echo "" >> final_report.txt

echo "SNP Distribution" \
>> final_report.txt

cat snp_statistics.txt >> final_report.txt

########################################################

########################################################
# STEP 15: DISPLAY RESULTS ON SCREEN
########################################################

echo ""
echo "========================================="
echo "           ANALYSIS SUMMARY"
echo "========================================="

echo ""
echo "===== FINAL REPORT ====="
cat final_report.txt

echo ""
echo "===== FIRST 10 STAT6 VARIANTS ====="
head STAT6_variant_table.tsv

echo ""
echo "===== FIRST 10 FCER1A VARIANTS ====="
head FCER1A_variant_table.tsv

echo ""
echo "===== SNP DISTRIBUTION ====="
cat snp_statistics.txt

echo ""
echo "========================================="
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================="

