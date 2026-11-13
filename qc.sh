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
Sortmerna_Root="$Outdir/4_sortmerna"
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

trim_adapters() {
    need "$cutadapt_txt" 1 || return 1

    header=$(head -n 1 "$cutadapt_txt")
    row=$(awk -F'\t' -v s="$sample_name" '$1 == s {print; exit}' "$cutadapt_txt")
    if [ -z "$row" ]; then
        log "ERROR: $sample_name not found in cutadapt.txt" && return 1
    fi

    IFS=$'\t' read -r -a cols <<< "$header"
    IFS=$'\t' read -r -a vals <<< "$row"
    r1_name=""
    r2_name=""
    cut_args=()
    for i in "${!cols[@]}"; do
        col="${cols[$i]}"
        val="${vals[$i]:-}"
        case "$col" in
            ""|sample_name) continue ;;
            r1)             r1_name="$val" ;;
            r2)             r2_name="$val" ;;
            -o|-p|-j)       log "WARN: column $col is set by the script, ignored" ;;
            -*)             [ -n "$val" ] && cut_args+=("$col" "$val") ;;
            *)              log "WARN: unknown column $col, ignored" ;;
        esac
    done

    if [ -z "$r1_name" ] || [ -z "$r2_name" ]; then
        log "ERROR: r1 and r2 columns must be set for $sample_name" && return 1
    fi

    r1="$Reads/$r1_name"
    r2="$Reads/$r2_name"
    need "$r1" 1 || return 1
    need "$r2" 1 || return 1

    skip_if_done "trim_adapters" "$trimmed_r1" "$trimmed_r2" && return 0
    mkdir -p "$sample_cutadapt"

    run cutadapt "${cut_args[@]}" \
        -j "$threads" \
        -o "$trimmed_r1" \
        -p "$trimmed_r2" \
        --report=full \
        "$r1" "$r2"
}

filter_rrna() {
    need "$trimmed_r1" 2 || return 1
    need "$trimmed_r2" 2 || return 1
    need "$Sortmerna_DB" || return 1
    skip_if_done "filter_rrna" "$unaligned" && return 0
    mkdir -p "$sample_sortmerna"
    run sortmerna \
        --ref "$Sortmerna_DB/silva-arc-16s-id95.fasta" \
        --ref "$Sortmerna_DB/silva-arc-23s-id98.fasta" \
        --ref "$Sortmerna_DB/silva-bac-16s-id90.fasta" \
        --ref "$Sortmerna_DB/silva-bac-23s-id98.fasta" \
        --ref "$Sortmerna_DB/silva-euk-18s-id95.fasta" \
        --ref "$Sortmerna_DB/silva-euk-28s-id98.fasta" \
        --ref "$Sortmerna_DB/rfam-5s-database-id98.fasta" \
        --ref "$Sortmerna_DB/rfam-5.8s-database-id98.fasta" \
        --workdir "$sample_sortmerna" \
        --aligned "$sample_sortmerna/${sample_name}_aligned" \
        --other "$sample_sortmerna/${sample_name}_unaligned" \
        --fastx \
        --out2 \
        --paired_out \
        --threads "$threads" \
        --reads "$trimmed_r1" \
        --reads "$trimmed_r2" || return 1
    cat \
        "$sample_sortmerna/${sample_name}_unaligned_fwd.fq" \
        "$sample_sortmerna/${sample_name}_unaligned_rev.fq" \
        > "$unaligned"
}

sample_done() {                             # usage: sample_done <sample_name>
    [ -e "$Cutadapt/$1/${1}_trimmed_R1.fastq" ] &&
    [ -e "$Cutadapt/$1/${1}_trimmed_R2.fastq" ] &&
    [ -e "$Sortmerna_Root/$1/${1}.fq" ]
}

trim_filter() {                                # asks before trimming, then loops over samples
    samples=$(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$cutadapt_txt")

    pending=0
    for s in $samples; do
        sample_done "$s" || pending=1
    done
    if [ "$pending" = 0 ] && [ "$force" != 1 ]; then
        log "Trimming and rRNA filtering already done for all samples"
        return 0
    fi

    printf 'Continue with trim_adapters and filter_rrna using the current cutadapt.txt? [y/n] '
    read -r answer || answer=n              # no input (non-interactive) counts as no
    case "$answer" in
        y|Y) ;;
        *)   log "Stopped before trimming. Edit cutadapt.txt, then run again."; return 0 ;;
    esac

    failed=""
    for sample_name in $samples; do
        log "=== Sample $sample_name ==="
        sample_cutadapt="$Cutadapt/$sample_name"
        sample_sortmerna="$Sortmerna_Root/$sample_name"
        trimmed_r1="$sample_cutadapt/${sample_name}_trimmed_R1.fastq"
        trimmed_r2="$sample_cutadapt/${sample_name}_trimmed_R2.fastq"
        unaligned="$sample_sortmerna/${sample_name}.fq"

        { trim_adapters && filter_rrna; } || failed="$failed $sample_name"
    done
    if [ -n "$failed" ]; then
        log "Failed samples:$failed"
        return 1
    fi
}

case "${1:-all}" in
    all)      run_stages 0 move_reads run_fastqc run_multiqc make_cutadapt_txt trim_filter ;;
    allforce) run_stages 1 move_reads run_fastqc run_multiqc make_cutadapt_txt trim_filter ;;
    *)        echo "Usage: $0 [all|allforce]"; exit 1 ;;
esac