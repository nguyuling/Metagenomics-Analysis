#!/bin/bash
set -e

# bash biobakery.sh all        run all unfinished steps
# bash biobakery.sh allforce   run all steps regardlessly

# helper function script reference
Scriptdir="$(cd "$(dirname "$0")" && pwd)"
source "$Scriptdir/lib.sh"

# global output directories and files
Reads="$Outdir/0_reads"
Cutadapt="$Outdir/3_cutadapt"
cutadapt_txt="$Outdir/3_cutadapt/cutadapt.txt"
Cross_Analysis="$Outdir/9_cross_analysis"

run_pair() {
    sample_name="$1"

    # per-sample output directories
    Cutadapt="$Outdir/3_cutadapt/$sample_name"
    Sortmerna="$Outdir/4_sortmerna/$sample_name"
    Metaphlan="$Outdir/5_metaphlan/$sample_name"
    Humann="$Outdir/6_humann/$sample_name"
    Genes_Mapping="$Outdir/7_genes_mapping/$sample_name"
    Pathways_Analysis="$Outdir/8_pathways_analysis/$sample_name"

    # per-sample output files
    trimmed_r1="$Cutadapt/${sample_name}_trimmed_R1.fastq"
    trimmed_r2="$Cutadapt/${sample_name}_trimmed_R2.fastq"
    unaligned="$Sortmerna/${sample_name}.fq"
    taxonomic_profile="$Metaphlan/${sample_name}_taxonomic_profile.tsv"
    genefamilies="$Humann/${sample_name}_genefamilies.tsv"
    genefamilies_renormed_uniref90="$Genes_Mapping/renormed/${sample_name}_renormed_uniref90.tsv"
    pathabundance="$Humann/${sample_name}_pathabundance.tsv"
    pathcoverage="$Humann/${sample_name}_pathcoverage.tsv"
    pathabundance_renormed="$Outdir/tmp/${sample_name}_pathabundance_renormed.tsv"
    pathabundance_unpacked="$Pathways_Analysis/${sample_name}_pathabundance_unpacked.tsv"
    pathcoverage_unpacked="$Pathways_Analysis/${sample_name}_pathcoverage_unpacked.tsv"

    log "=== Sample $sample_name ==="
    for stage in trim_adapters filter_rrna taxonomic_profiling \
                 functional_profiling genes_mapping unpack_pathways; do
        "$stage" || { log "Sample $sample_name stopped at $stage"; return 1; }
    done
}

run_all_pairs() {                           # usage: run_all_pairs <force 0|1>
    force=$1
    before=${ran_count:-0}
    failed=""

    for s in $(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$cutadapt_txt"); do
        run_pair "$s" || failed="$failed $s"
    done

    if [ -n "$failed" ]; then
        log "Step 8 skipped: failed samples:$failed"
        return 1
    fi

    if [ "${ran_count:-0}" -gt "$before" ] || [ "$force" = 1 ] || [ ! -e "$Cross_Analysis/genefamilies_renormed.tsv" ]; then
        cross_analysis || { log "Cross-analysis failed"; return 1; }
    else
        log "Step 8 skipped: no upstream changes"
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
    mkdir -p "$Cutadapt"

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
    mkdir -p "$Sortmerna"
    run sortmerna \
        --ref "$Sortmerna_DB/silva-arc-16s-id95.fasta" \
        --ref "$Sortmerna_DB/silva-arc-23s-id98.fasta" \
        --ref "$Sortmerna_DB/silva-bac-16s-id90.fasta" \
        --ref "$Sortmerna_DB/silva-bac-23s-id98.fasta" \
        --ref "$Sortmerna_DB/silva-euk-18s-id95.fasta" \
        --ref "$Sortmerna_DB/silva-euk-28s-id98.fasta" \
        --ref "$Sortmerna_DB/rfam-5s-database-id98.fasta" \
        --ref "$Sortmerna_DB/rfam-5.8s-database-id98.fasta" \
        --workdir "$Sortmerna" \
        --aligned "$Sortmerna/${sample_name}_aligned" \
        --other "$Sortmerna/${sample_name}_unaligned" \
        --fastx \
        --out2 \
        --paired_out \
        --threads "$threads" \
        --reads "$trimmed_r1" \
        --reads "$trimmed_r2"  || return 1
        # SET paired_in / paired_out (aligned if either read align / only if both reads align)
    cat \
        "$Sortmerna/${sample_name}_unaligned_fwd.fq" \
        "$Sortmerna/${sample_name}_unaligned_rev.fq" \
        > "$unaligned"
}

taxonomic_profiling() {
    need "$unaligned" 3 || return 1
    need "$Metaphlan_DB" || return 1
    skip_if_done "taxonomic_profiling" "$taxonomic_profile" "$Metaphlan/$sample_name.krona.html" && return 0
    mkdir -p "$Metaphlan"
    run metaphlan \
        --input_type fastq \
        --bowtie2db "$Metaphlan_DB" \
        -x mpa_vJun23_CHOCOPhlAnSGB_202307 \
        --bowtie2out "$Metaphlan/${sample_name}_metagenome_vJun23.bowtie2.bz2" \
        -o "$taxonomic_profile" \
        --nproc "$threads" \
        --offline \
        "$unaligned"  || return 1
    run metaphlan2krona.py \
        -p "$taxonomic_profile" \
        -k "$Metaphlan/${sample_name}_krona.txt"  || return 1
    run ktImportText \
        -o "$Metaphlan/$sample_name.krona.html" \
        -n 'all' \
        "$Metaphlan/${sample_name}_krona.txt"
        # SET n (highest taxa level name)
}

functional_profiling() {
    need "$unaligned" 3 || return 1
    need "$taxonomic_profile" 4 || return 1
    need "$Humann_DB/chocophlan" || return 1
    need "$Humann_DB/uniref" || return 1
    skip_if_done "functional_profiling" "$genefamilies" "$pathabundance" "$pathcoverage" && return 0
    mkdir -p "$Humann"
    run humann \
        -i "$unaligned" \
        -o "$Humann" \
        --taxonomic-profile "$taxonomic_profile" \
        --nucleotide-database "$Humann_DB/chocophlan" \
        --protein-database "$Humann_DB/uniref" \
        --threads "$threads"
}

sort_table() {
    file="$1"
    header=$(head -n 1 "$file")
    tail -n +2 "$file" | sort -t "$(printf '\t')" -k2,2gr > "$file.tmp"
    { printf '%s\n' "$header"; cat "$file.tmp"; } > "$file"
    rm -f "$file.tmp"
}

genes_mapping() {
    need "$genefamilies" 5 || return 1
    need "$Humann_DB/utility_mapping" || return 1
    path="$Genes_Mapping/${sample_name}"
    skip_if_done "genes_mapping" "${path}_uniref90.tsv" "${path}_uniref50.tsv" "${path}_eggnog.tsv" \
        "${path}_go.tsv" "${path}_ko.tsv" "${path}_level4ec.tsv" "${path}_pfam.tsv" && return 0
    mkdir -p "$Genes_Mapping/regrouped" "$Genes_Mapping/renamed" "$Genes_Mapping/renormed"

    while read -r category groupmap namemap; do
        regrouped="$Genes_Mapping/regrouped/${sample_name}_regrouped_${category}.tsv"
        renamed="$Genes_Mapping/renamed/${sample_name}_renamed_${category}.tsv"
        renormed="$Genes_Mapping/renormed/${sample_name}_renormed_${category}.tsv"
        combined="$Genes_Mapping/${sample_name}_${category}.tsv"

        if [ "$groupmap" = "none" ]; then
            cp "$genefamilies" "$regrouped"  || return 1
        else
            run humann_regroup_table \
                -i "$genefamilies" \
                -c "$Humann_DB/utility_mapping/$groupmap" \
                -o "$regrouped" || return 1
        fi
        run sort_table "$regrouped" || return 1
        run humann_rename_table \
            -i "$regrouped" \
            -c "$Humann_DB/utility_mapping/$namemap" \
            -o "$renamed" || return 1
        run humann_renorm_table \
            -i "$regrouped" \
            -u relab \
            -o "$renormed" || return 1
        paste "$renamed" \
            <(awk -F'\t' 'BEGIN{OFS="\t"} {print $NF}' "$renormed") > "$combined" || return 1

    done <<EOF
uniref90    none                            map_uniref90_name.txt.bz2
uniref50    map_uniref50_uniref90.txt.gz    map_uniref50_name.txt.bz2
eggnog      map_eggnog_uniref90.txt.gz      map_eggnog_name.txt.gz
go          map_go_uniref90.txt.gz          map_go_name.txt.gz
ko          map_ko_uniref90.txt.gz          map_ko_name.txt.gz
level4ec    map_level4ec_uniref90.txt.gz    map_ec_name.txt.gz
pfam        map_pfam_uniref90.txt.gz        map_pfam_name.txt.gz
EOF

    rm -f "$Outdir/$sample_name.xlsx"          # start a fresh workbook for this sample
    run python "$Scriptdir/make_excel.py" "$sample_name" "$Outdir" "$Genes_Mapping" || return 1
}

unpack_pathways() {
    need "$genefamilies_renormed_uniref90" 6 || return 1
    need "$pathabundance" 5 || return 1
    need "$pathcoverage" 5 || return 1
    need "$Humann_DB/pathways" || return 1
    skip_if_done "unpack_pathways" "$pathabundance_unpacked" "$pathcoverage_unpacked" && return 0
    mkdir -p "$Pathways_Analysis" "$Outdir/tmp"
    rm -f "$pathabundance_renormed"                 # clear any leftover from an interrupted run

    run sort_table "$pathabundance" || return 1
    run humann_renorm_table \
        -i "$pathabundance" \
        -u relab \
        -o "$pathabundance_renormed" || return 1
    run humann_unpack_pathways \
        --input-genes "$genefamilies_renormed_uniref90" \
        --input-pathways "$pathabundance_renormed" \
        --gene-mapping "$Humann_DB/pathways/metacyc_reactions_level4ec_only.uniref.bz2" \
        --pathway-mapping "$Humann_DB/pathways/metacyc_pathways_structured_filtered_v24" \
        --output "$pathabundance_unpacked" || return 1
    rm -f "$pathabundance_renormed"                 # temporary file, no longer needed

    run sort_table "$pathcoverage" || return 1
    run humann_unpack_pathways \
        --input-genes "$genefamilies_renormed_uniref90" \
        --input-pathways "$pathcoverage" \
        --gene-mapping "$Humann_DB/pathways/metacyc_reactions_level4ec_only.uniref.bz2" \
        --pathway-mapping "$Humann_DB/pathways/metacyc_pathways_structured_filtered_v24" \
        --output "$pathcoverage_unpacked" || return 1
    run python "$Scriptdir/make_excel.py" "$sample_name" "$Outdir" "$Pathways_Analysis" || return 1
}

cross_analysis() {
    skip_if_done "cross_analysis" "$Cross_Analysis/genefamilies_renormed.tsv" \
        "$Cross_Analysis/pathabundance.tsv" "$Cross_Analysis/pathcoverage.tsv" \
        "$Outdir/cross_analysis.xlsx" && return 0
    mkdir -p "$Cross_Analysis"
    samples=$(awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "$cutadapt_txt")
    stage="$Cross_Analysis/.stage"
    rm -rf "$stage"
    mkdir -p "$stage/genefamilies" "$stage/pathabundance" "$stage/pathcoverage"
    for s in $samples; do
        for t in genefamilies pathabundance pathcoverage; do
            src="$Outdir/6_humann/$s/${s}_${t}.tsv"
            [ -e "$src" ] || { log "ERROR: missing $src"; rm -rf "$stage"; return 1; }
            ln -sf "$src" "$stage/$t/${s}_${t}.tsv"
        done
    done

    run humann_join_tables -i "$stage/genefamilies" \
        -o "$Cross_Analysis/genefamilies.tsv" || { rm -rf "$stage"; return 1; }
    run humann_rename_table \
        -i "$Cross_Analysis/genefamilies.tsv" \
        -c "$Humann_DB/utility_mapping/map_uniref90_name.txt.bz2" \
        -o "$Cross_Analysis/genefamilies.tmp" \
    && mv "$Cross_Analysis/genefamilies.tmp" "$Cross_Analysis/genefamilies.tsv" \
    || { rm -rf "$stage"; return 1; }
    run humann_renorm_table \
        -i "$Cross_Analysis/genefamilies.tsv" -u relab \
        -o "$Cross_Analysis/genefamilies_renormed.tsv" || { rm -rf "$stage"; return 1; }
    run humann_join_tables -i "$stage/pathabundance" \
        -o "$Cross_Analysis/pathabundance.tsv" || { rm -rf "$stage"; return 1; }
    run humann_join_tables -i "$stage/pathcoverage" \
        -o "$Cross_Analysis/pathcoverage.tsv" || { rm -rf "$stage"; return 1; }
    rm -rf "$stage"
    run humann_renorm_table \
        -i "$Cross_Analysis/pathabundance.tsv" -u relab \
        -o "$Cross_Analysis/pathabundance_renormed.tsv" || return 1
    run python "$Scriptdir/make_excel.py" "cross_analysis" "$Outdir" "$Cross_Analysis" || return 1
}

# entry point
case "${1:-}" in
    all)        run_all_pairs 0 ;;
    allforce)   run_all_pairs 1 ;;
    9)          force=1; generate_barplot "$2" "$3" || exit 1 ;;
    *)          echo "Usage: bash $0 [all   |   allforce    |   9 <table_path> <gene_id>]"; exit 1 ;;
esac