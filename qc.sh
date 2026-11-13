#!/bin/bash
set -e

# bash qc.sh all        run all unfinished steps
# bash qc.sh allforce   run all steps regardlessly

# helper function script reference
Scriptdir="$(cd "$(dirname "$0")" && pwd)"
source "$Scriptdir/lib.sh"

# global output directories and files
Reads="$Outdir/0_reads"
Fastqc="$Outdir/1_fastqc"
Multiqc="$Outdir/2_multiqc"
Cutadapt="$Outdir/3_cutadapt"
cutadapt_txt="$Cutadapt/cutadapt.txt"

move_reads() {
    mkdir -p "$Reads"
    if ls "$Outdir"/*.fastq* >/dev/null 2>&1; then  # if read files exist in outdir
        mv "$Outdir"/*.fastq* "$Reads" || return 1  # move them to reads folder                    
        log "Moved reads to $Reads"
    else
        log "No .fastq* files in $Outdir, skipping move."   # handle subsequent runs
    fi
}

run_fastqc() {
    need "$Reads" 1 || return 1
    run mkdir -p "$Fastqc" || return 1
    for f in "$Reads"/*.fastq*; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        name=${base%.fastq.gz}          # strip .fastq.gz first
        name=${name%.fastq}             # then plain .fastq
        skip_if_done "fastqc $base" "$Fastqc/${name}_fastqc.html" && continue
        run fastqc -o "$Fastqc" -t "$threads" "$f" || return 1
    done
}

run_multiqc() {
    need "$Fastqc" 1 || return 1
    skip_if_done "multiqc" "$Multiqc/multiqc_report.pdf" && return 0
    mkdir -p "$Multiqc"
    run multiqc \
        -n multiqc_report \
        -o "$Multiqc" \
        --pdf "$Fastqc" || return 1
}

make_cutadapt_txt() {
    need "$Reads" 1 || return 1
    skip_if_done "make_cutadapt_txt" "$cutadapt_txt" && return 0
    mkdir -p "$(dirname "$cutadapt_txt")"
    mkdir -p "$Cutadapt"

    # mapping table: column_header     value
        mapping='sample_name    "$sample_name"
    r1    "$r1"
    r2    "$r2"
    -a    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    -A    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    -q    20
    -m    20'

    header=$(printf '%s\n' "$mapping" | while read -r header value; do  # read header and value columns from mapping
        printf '%s\t' "$header"                         # append headers separated by tab
    done)
    printf '%s\n' "${header%\t}" > "$cutadapt_txt"      # write header to cutadapt.txt

    paired=0                                            # init counters
    unpaired=0

    for f in "$Reads"/*_R1*; do                         # for each read file containing _R1
        [ -f "$f" ] || continue                         # check if it is file
        r1=$(basename "$f")                             # store r1 file name
        sample_name="${r1%%_R1*}"                       # sample name = text before _R1
        r2=$(echo "$r1" | sed 's/_R1/_R2/')             # generate expected r2 filename
        if [ ! -f "$Reads/$r2" ]; then                  # if expected r2 not found
            unpaired=$((unpaired + 1))                  # increase unpaired counter
            continue                                    # continue to find pair reads
        fi
        row=$(printf '%s\n' "$mapping" | while read -r header value; do # read header and value column from mapping
            eval "v=$value"                             # expand value (to handle reference variables)
            printf '%s\t' "$v"                          # append values separated by tab
        done)
        printf '%s\n' "${row%$\t}" >> "$cutadapt_txt"   # write values to cutadapt.txt
        paired=$((paired + 1))                          # increase paired counter
    done

    for f in "$Reads"/*_R2*; do                         # for each read file containing _R2
        [ -f "$f" ] || continue                         # check if it is file
        r2=$(basename "$f")                             # store r2 file name
        r1=$(echo "$r2" | sed 's/_R2/_R1/')             # generate expected r1 filename
        [ -f "$Reads/$r1" ] || unpaired=$((unpaired + 1))   # if expected r1 not found then increase unpaired counter
    done

    log "$paired pair(s) written to cutadapt.txt, edit the flag (-a, -A, -q, -m) if needed."
    if [ "$unpaired" -gt 0 ]; then
        log "WARNING: $unpaired file(s) have no matching pair and are NOT included in cutadapt.txt."
        log "Please add their sample_name, r1 and r2 to cutadapt.txt manually."
    fi
}

case "${1:-all}" in
    all)      run_stages 0 move_reads run_fastqc run_multiqc make_cutadapt_txt ;;
    allforce) run_stages 1 move_reads run_fastqc run_multiqc make_cutadapt_txt ;;
    *)        echo "Usage: $0 [all|allforce]"; exit 1 ;;
esac