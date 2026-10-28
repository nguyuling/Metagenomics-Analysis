#!/bin/sh
set -e

# RENAME
Outdir=/Volumes/T7/metagenomics_project
Refdir=/Volumes/T7/databases/sortmerna
r1=T1A_R1.fastq
r2=T1A_R2.fastq
threads=8

# directories and files
Reads="$Outdir/1_reads"
Fastqc="$Outdir/2_fastqc"
Multiqc="$Outdir/3_multiqc"
Trimmed="$Outdir/4_trimmed"
Sortmerna="$Outdir/5_sortmerna"

r1_trimmed="${r1%.fastq*}_trimmed.${r1##*.}"
r2_trimmed="${r2%.fastq*}_trimmed.${r2##*.}"