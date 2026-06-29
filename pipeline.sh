#!/bin/bash

########################################################
# ADVANCED FOOD ALLERGY VARIANT ANALYSIS PIPELINE
########################################################

echo "========================================="
echo " FOOD ALLERGY NGS PIPELINE STARTED"
echo "========================================="

REF=STAT6_FCER1A.fasta
THREADS=8

########################################################
# MULTI-SAMPLE SUPPORT
########################################################

for SAMPLE in $(cat samples.txt)
do

echo ""
echo "========================================="
echo "Processing Sample: $SAMPLE"
echo "========================================="

########################################################
# STEP 1: DOWNLOAD DATA
########################################################

if [ ! -f "${SAMPLE}_1.fastq" ]; then

    echo "Downloading SRA sample..."

    prefetch $SAMPLE

    fasterq-dump $SAMPLE \
    -e $THREADS

fi

########################################################
# STEP 2: FASTQC
########################################################

fastqc ${SAMPLE}_1.fastq
fastqc ${SAMPLE}_2.fastq

########################################################
# STEP 3: TRIMMING
########################################################

fastp \
-i ${SAMPLE}_1.fastq \
-I ${SAMPLE}_2.fastq \
-o ${SAMPLE}_1_clean.fastq \
-O ${SAMPLE}_2_clean.fastq \
--cut_right \
--cut_right_mean_quality 20 \
--length_required 30 \
-w $THREADS \
-h ${SAMPLE}_fastp.html \
-j ${SAMPLE}_fastp.json

########################################################
# STEP 4: FASTQC AFTER CLEANING
########################################################

fastqc ${SAMPLE}_1_clean.fastq
fastqc ${SAMPLE}_2_clean.fastq

########################################################
# STEP 5: INDEX REFERENCE
########################################################

if [ ! -f "${REF}.bwt" ]; then

    bwa index $REF

    samtools faidx $REF

fi

########################################################
# STEP 6: ALIGNMENT
########################################################

bwa mem \
-t $THREADS \
$REF \
${SAMPLE}_1_clean.fastq \
${SAMPLE}_2_clean.fastq | \
samtools sort \
-@ $THREADS \
-o ${SAMPLE}_sorted.bam

########################################################
# STEP 7: BAM INDEX
########################################################

samtools index ${SAMPLE}_sorted.bam

########################################################
# STEP 8: VARIANT CALLING
########################################################

bcftools mpileup \
-f $REF \
${SAMPLE}_sorted.bam | \
bcftools call \
-mv \
-Ov \
-o ${SAMPLE}_variants.vcf

########################################################
# STEP 9: VARIANT ANNOTATION (VEP)
########################################################

vep \
-i ${SAMPLE}_variants.vcf \
-o ${SAMPLE}_annotated.txt \
--cache \
--assembly GRCh38 \
--sift b \
--polyphen b

########################################################
# STEP 10: GENE-WISE EXTRACTION
########################################################

grep "NC_000012" \
${SAMPLE}_variants.vcf \
> ${SAMPLE}_STAT6_variants.txt

grep "NC_000001" \
${SAMPLE}_variants.vcf \
> ${SAMPLE}_FCER1A_variants.txt

STAT6_COUNT=$(grep -v "^#" ${SAMPLE}_STAT6_variants.txt | wc -l)

FCER1A_COUNT=$(grep -v "^#" ${SAMPLE}_FCER1A_variants.txt | wc -l)

TOTAL=$(grep -v "^#" ${SAMPLE}_variants.vcf | wc -l)

########################################################
# STEP 11: VARIANT TABLES
########################################################

echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \
> ${SAMPLE}_STAT6_table.tsv

grep -v "^#" \
${SAMPLE}_STAT6_variants.txt | \
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
>> ${SAMPLE}_STAT6_table.tsv

echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \
> ${SAMPLE}_FCER1A_table.tsv

grep -v "^#" \
${SAMPLE}_FCER1A_variants.txt | \
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
>> ${SAMPLE}_FCER1A_table.tsv

########################################################
# STEP 12: SNP STATISTICS
########################################################

grep -v "^#" ${SAMPLE}_variants.vcf | \
awk '{print $4">"$5}' | \
sort | uniq -c \
> ${SAMPLE}_snp_statistics.txt

########################################################
# STEP 13: GENERATE PLOT
########################################################

echo "STAT6 $STAT6_COUNT" > counts.txt
echo "FCER1A $FCER1A_COUNT" >> counts.txt

gnuplot << EOF

set terminal png
set output "${SAMPLE}_variant_counts.png"

set style data histograms
set style fill solid

set xlabel "Genes"
set ylabel "Variant Count"

plot 'counts.txt' using 2:xtic(1) title 'Variants'

EOF

########################################################
# STEP 14: FINAL REPORT
########################################################

echo "Food Allergy Variant Analysis Report" \
> ${SAMPLE}_final_report.txt

echo "==================================" \
>> ${SAMPLE}_final_report.txt

echo "" >> ${SAMPLE}_final_report.txt

echo "Sample ID: $SAMPLE" \
>> ${SAMPLE}_final_report.txt

echo "Total Variants: $TOTAL" \
>> ${SAMPLE}_final_report.txt

echo "STAT6 Variants: $STAT6_COUNT" \
>> ${SAMPLE}_final_report.txt

echo "FCER1A Variants: $FCER1A_COUNT" \
>> ${SAMPLE}_final_report.txt

echo "" >> ${SAMPLE}_final_report.txt

echo "SNP Distribution" \
>> ${SAMPLE}_final_report.txt

cat ${SAMPLE}_snp_statistics.txt \
>> ${SAMPLE}_final_report.txt

########################################################
# STEP 15: DISPLAY RESULTS
########################################################

echo ""
echo "===== REPORT ====="
cat ${SAMPLE}_final_report.txt

echo ""
echo "===== FIRST 10 STAT6 VARIANTS ====="
head ${SAMPLE}_STAT6_table.tsv

echo ""
echo "===== FIRST 10 FCER1A VARIANTS ====="
head ${SAMPLE}_FCER1A_table.tsv

########################################################
# STEP 16: CLEANUP
########################################################

rm -f *.fastq
rm -f *.sra

echo "Temporary files removed."

done

########################################################
# STEP 17: MULTIQC REPORT
########################################################

multiqc . -o multiqc_output

echo ""
echo "========================================="
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================="
