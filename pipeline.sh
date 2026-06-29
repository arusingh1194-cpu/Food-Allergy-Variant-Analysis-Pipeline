#!/bin/bash

########################################################

# FOOD ALLERGY VARIANT ANALYSIS PIPELINE

# Optimized for WSL/Laptop

########################################################

echo "========================================="
echo " FOOD ALLERGY NGS PIPELINE STARTED"
echo "========================================="

REF=STAT6_FCER1A.fasta
THREADS=4

########################################################

# REFERENCE INDEX

########################################################

if [ ! -f "${REF}.bwt" ]; then

```
echo "Indexing reference..."

bwa index $REF
samtools faidx $REF
```

fi

########################################################

# PROCESS ALL SAMPLES

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

```
echo "Downloading sample..."

prefetch $SAMPLE

fasterq-dump $SAMPLE \
-X 100000 \
-e $THREADS \
-p
```

fi

########################################################

# STEP 2: FASTQC

########################################################

if [ ! -f "${SAMPLE}_1_fastqc.html" ]; then

```
fastqc ${SAMPLE}_1.fastq
fastqc ${SAMPLE}_2.fastq
```

fi

########################################################

# STEP 3: TRIMMING

########################################################

if [ ! -f "${SAMPLE}_1_clean.fastq" ]; then

```
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
```

fi

########################################################

# STEP 4: ALIGNMENT

########################################################

if [ ! -f "${SAMPLE}_sorted.bam" ]; then

```
bwa mem \
-t $THREADS \
$REF \
${SAMPLE}_1_clean.fastq \
${SAMPLE}_2_clean.fastq | \
samtools sort \
-@ $THREADS \
-o ${SAMPLE}_sorted.bam
```

fi

########################################################

# STEP 5: BAM INDEX

########################################################

if [ ! -f "${SAMPLE}_sorted.bam.bai" ]; then

```
samtools index ${SAMPLE}_sorted.bam
```

fi

########################################################

# STEP 6: VARIANT CALLING

########################################################

if [ ! -f "${SAMPLE}_variants.vcf" ]; then

```
bcftools mpileup \
-f $REF \
${SAMPLE}_sorted.bam | \
bcftools call \
-mv \
-Ov \
-o ${SAMPLE}_variants.vcf
```

fi

########################################################

# STEP 7: GENE EXTRACTION

########################################################

grep "NC_000012" 
${SAMPLE}_variants.vcf \

> ${SAMPLE}_STAT6.txt

grep "NC_000001" 
${SAMPLE}_variants.vcf \

> ${SAMPLE}_FCER1A.txt

STAT6_COUNT=$(grep -v "^#" ${SAMPLE}_STAT6.txt | wc -l)

FCER1A_COUNT=$(grep -v "^#" ${SAMPLE}_FCER1A.txt | wc -l)

TOTAL=$(grep -v "^#" ${SAMPLE}_variants.vcf | wc -l)

########################################################

# STEP 8: VARIANT TABLES

########################################################

echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \

> ${SAMPLE}_STAT6_table.tsv

grep -v "^#" ${SAMPLE}_STAT6.txt | 
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \

> > ${SAMPLE}_STAT6_table.tsv

echo -e "Chromosome\tPosition\tREF\tALT\tQUAL" \

> ${SAMPLE}_FCER1A_table.tsv

grep -v "^#" ${SAMPLE}_FCER1A.txt | 
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \

> > ${SAMPLE}_FCER1A_table.tsv

########################################################

# STEP 9: SNP STATISTICS

########################################################

grep -v "^#" ${SAMPLE}_variants.vcf | 
awk '{print $4">"$5}' | 
sort | uniq -c \

> ${SAMPLE}_snp_statistics.txt

########################################################

# STEP 10: PLOT GENERATION

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

# STEP 11: FINAL REPORT

########################################################

echo "Food Allergy Variant Analysis Report" \

> ${SAMPLE}_final_report.txt

echo "===================================" \

> > ${SAMPLE}_final_report.txt

echo "" >> ${SAMPLE}_final_report.txt

echo "Sample ID: $SAMPLE" \

> > ${SAMPLE}_final_report.txt

echo "Total Variants: $TOTAL" \

> > ${SAMPLE}_final_report.txt

echo "STAT6 Variants: $STAT6_COUNT" \

> > ${SAMPLE}_final_report.txt

echo "FCER1A Variants: $FCER1A_COUNT" \

> > ${SAMPLE}_final_report.txt

echo "" >> ${SAMPLE}_final_report.txt

echo "SNP Distribution" \

> > ${SAMPLE}_final_report.txt

cat ${SAMPLE}_snp_statistics.txt \

> > ${SAMPLE}_final_report.txt

########################################################

# STEP 12: DISPLAY RESULTS

########################################################

echo ""
echo "===== REPORT ====="
cat ${SAMPLE}_final_report.txt

echo ""
echo "===== STAT6 VARIANTS ====="
head ${SAMPLE}_STAT6_table.tsv

echo ""
echo "===== FCER1A VARIANTS ====="
head ${SAMPLE}_FCER1A_table.tsv

########################################################

# STEP 13: CLEANUP

########################################################

rm -f ${SAMPLE}_1.fastq
rm -f ${SAMPLE}_2.fastq
rm -rf ${SAMPLE}

echo "Temporary files removed."

done

########################################################

# STEP 14: MULTIQC

########################################################

multiqc . -o multiqc_output

echo ""
echo "========================================="
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================="

echo "========================================="
