#!/bin/sh
set -e

mkdir -p "$Reads"
mv "$Outdir"/*.fastq* "$Reads"/

mkdir -p "$Fastqc"
fastqc \
    -o "$Fastqc" \
    -t "$threads" \
    "$Reads"/*.fastq*

mkdir -p "$Multiqc"
multiqc \
    -n multiqc_report \
    -o "$Multiqc" \
    --pdf \
    "$Fastqc"

mkdir -p "$Trimmed"
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20 \
    -m 20 \
    -j 8 \
    -o "$Trimmed/$r1_trimmed" \
    -p "$Trimmed/$r2_trimmed" \
    --report=full \
    "$Reads/$r1" "$Reads/$r2" > "$Trimmed/cutadapt_report.txt"
# SET a (3' end adapter R1)
# SET A (3' end adapter R2)
# SET q (quality cutoff)
# SET m (minimum length)
# SET j (cores)

mkdir "$Sortmerna"
sortmerna \
    --ref "$Refdir/silva-arc-16s-id95.fasta" \
    --ref "$Refdir/silva-arc-23s-id98.fasta" \
    --ref "$Refdir/silva-bac-16s-id90.fasta" \
    --ref "$Refdir/silva-bac-23s-id98.fasta" \
    --ref "$Refdir/silva-euk-18s-id95.fasta" \
    --ref "$Refdir/silva-euk-28s-id98.fasta" \
    --ref "$Refdir/rfam-5s-database-id98.fasta" \
    --ref "$Refdir/rfam-5.8s-database-id98.fasta" \
    --workdir "$Sortmerna" \
    --aligned "$Sortmerna/aligned" \
    --other "$Sortmerna/unaligned" \
    --fastx \
    --out2 \
    --paired_out \
    --threads "$threads" \
    --reads "$Trimmed/$r1_trimmed" \
    --reads "$Trimmed/$r2_trimmed"
# SET / ADD ref <Refdir/*.fasta>
# SET task (0 - align only)
# SET paired_in (aligned if either read align) / paired_out (aligned if both read align)