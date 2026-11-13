#!/bin/bash
set -e

# SET output and database directories
Outdir=/Volumes/T7/metagenomics_project             # output directory
Sortmerna_DB=/Volumes/T7/databases/sortmerna        # rrna sequences database
Metaphlan_DB=/Volumes/T7/conda_envs/biobakery/lib/python3.10/site-packages/metaphlan/metaphlan_databases # metaphlan database
Humann_DB=/Volumes/T7/conda_envs/biobakery/lib/python3.10/site-packages/humann/data
Log="$Outdir/logs"                                  # log directories
threads=8                                           # number of core

mkdir -p "$Log"
log_file="$Log/pipeline_$(date '+%Y%m%d_%H%M%S').log"

log() {                                     # log datetime and message to terminal and log file
    msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$log_file"
}

run() {                                     # usage: run <command> <parameter>
    log "START: $*"                         # log the starting step
    start=$(date +%s)                       # start time
    if "$@" >> "$log_file" 2>&1; then       # run the command, write stdout / stderr to log file
        ran_count=$(( ${ran_count:-0} + 1 ))    
        log "OK ($(( $(date +%s) - start ))s)"  # exit 0 : log run time (current time - start time)
    else
        rc=$?                               # exit != 0 : create variable of the latest exit code
        log "FAILED (exit $rc): $*"         # log the exit code and failed command
        return "$rc"
    fi
}

need() {                                    # usage: need <file or dir> <step number that creates it>
    if [ ! -e "$1" ]; then                  # if the prerequisite file is absent
        log "ERROR: missing $1. Run step $2 first."  && return 1   # log reminder message about running the upstream step first
    fi
}

skip_if_done() {                            # usage: skip_if_done <step name> <output>...
    name=$1; shift                          # remove first arguments (step name), leaving only the output files
    [ "${force:-0}" = 1 ] && { log "RUN $name: force=1"; return 1; }
    for f in "$@"; do
        if [ ! -e "$f" ]; then
            log "RUN $name: missing $f"
            return 1
        fi
    done
    log "SKIPPED $name: outputs already exist"
    return 0
}

run_stages() {                              # usage: run_stages <force 0|1> <function>...
    force=$1; shift
    start=$(date +%s)
    for stage in "$@"; do
        if ! "$stage"; then
            log "STOPPED at $stage. Fix the error and run again."
            return 1
        fi
    done
    total=$(( $(date +%s) - start ))
    log "DONE. Total time: $(( total / 60 )) min $(( total % 60 )) sec"
}