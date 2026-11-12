# shellcheck shell=bash
# Sourced by every bash cell of benchmarks.qmd: settings, then `bench`.
#
#   bench SECTION DATASET TOOL TASK OUTPUT -- COMMAND [ARG...]
#
# runs COMMAND as a fresh process, `reps DATASET` times, each under GNU time
# and `timeout`, after deleting OUTPUT (the file it writes, '-' for none) and
# its sidecars; and appends one row per run to $RESULTS/timings.tsv, with its
# status: ok, timeout (stopped at $TIMEOUT s), oom (killed by the kernel for
# want of memory) or error. A run that was already timed, and whose output is
# still there, is not run again (BENCH_REDO=1 forces it), so a page can be
# re-rendered without re-timing.

set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${BENCH_WORK:-$HERE/_work}         # inputs, maps, logs (large)
RESULTS=${BENCH_RESULTS:-$HERE/_results} # tables and figures (small)
DATA=$WORK/data       # one directory per dataset: simulate.py's output
MAPS=$WORK/maps       # the maps later sections read
QUERIES=$WORK/queries # what the queries returned, to compare the two tools
SCRATCH=$WORK/scratch # outputs that are only timed, deleted after their runs
LOGS=$WORK/logs
TIMINGS=$RESULTS/timings.tsv

DATASETS=${BENCH_DATASETS:-S M L XL}
BINSIZE=${BENCH_BINSIZE:-1000}
TIMEOUT=${BENCH_TIMEOUT:-1800}         # seconds, per run
COOLER_VERSION=0.12.1                  # the cooler every timing is against; CI installs it
ZOOMS=${BENCH_RESOLUTIONS:-$BINSIZE,$((2 * BINSIZE)),$((5 * BINSIZE)),$((10 * BINSIZE))}

mkdir -p "$DATA" "$MAPS" "$QUERIES" "$SCRATCH" "$LOGS" "$RESULTS"

# Runs per command: more where a run is cheap. BENCH_REPS sets them all.
reps() {
    if [[ -n ${BENCH_REPS:-} ]]; then echo "$BENCH_REPS"; return; fi
    case $1 in S) echo 5 ;; M) echo 3 ;; L) echo 2 ;; *) echo 1 ;; esac
}

COLUMNS_TSV="section	dataset	tool	task	rep	status	wall_s	cpu_s	user_s	sys_s	peak_mb	output_mb	started	metacooler	cooler	cpu	cores	memory_gb	command"
[[ -s $TIMINGS ]] || echo "$COLUMNS_TSV" >"$TIMINGS"

# Versions and machine, once per work directory: every row carries them.
if [[ ! -s $WORK/.context ]]; then
    {
        python -c 'import metacooler, cooler; print(metacooler.__version__, cooler.__version__, sep="\t")' \
            2>/dev/null || printf 'NA\tNA\n'
        lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | grep . || echo NA
        nproc
        awk '/MemTotal/ {printf "%.1f\n", $2 / 1048576}' /proc/meminfo 2>/dev/null || echo NA
    } | head -n 4 | paste -sd '\t' >"$WORK/.context"
fi
CONTEXT=$(cat "$WORK/.context")

# Rows already in the table for this command, and the status of the last one.
_done() { awk -F'\t' -v k="$1" '$1"/"$2"/"$3"/"$4 == k { n++; s = $6 } END { print n + 0, s }' "$TIMINGS"; }

bench() {
    local section=$1 dataset=$2 tool=$3 task=$4 output=$5
    shift 5
    [[ ${1-} == -- ]] && shift
    local key="$section/$dataset/$tool/$task" n_reps n_done last first
    n_reps=$(reps "$dataset")
    if [[ ${BENCH_REDO:-0} == 1 ]]; then
        awk -F'\t' -v k="$key" '$1"/"$2"/"$3"/"$4 != k' "$TIMINGS" >"$TIMINGS.tmp"
        mv "$TIMINGS.tmp" "$TIMINGS"
    fi
    read -r n_done last <<<"$(_done "$key")"
    first=$((n_done + 1))
    if [[ $last == timeout || $last == oom ]]; then # final: it would only fail again
        printf '%-44s did not finish before (%s): not run again\n' "$key" "$last"
        return 0
    fi
    if ((n_done >= n_reps)); then
        if [[ $output == - || $output == "$SCRATCH"/* || -e $output || $last != ok ]]; then
            printf '%-44s already timed (%d runs, last %s)\n' "$key" "$n_done" "$last"
            return 0
        fi
        n_reps=$first # its output is gone: one more run to make it again
    fi
    local log="$LOGS/${key//\//_}" command="$*" rep status signal wall user sys rss bytes started
    command=${command//$'\t'/ }
    command=${command//$'\n'/ }
    for ((rep = first; rep <= n_reps; rep++)); do
        # The output, and the P(s) table `create pairs` writes beside it.
        [[ $output != - ]] && rm -rf "$output" "$output".* "${output%.*}.ps.tsv"
        started=$(date -u +%Y-%m-%dT%H:%M:%S)
        /usr/bin/time -f '%e %U %S %M %x' -o "$log.time" \
            timeout --kill-after=60 "$TIMEOUT" "$@" >"$log.out" 2>"$log.err"
        read -r wall user sys rss status <<<"$(tail -n 1 "$log.time")"
        # A killed command has no exit status (GNU time's %x reads 0): take
        # the signal from the line time writes about it, the shell's way.
        signal=$(sed -n 's/^Command terminated by signal \([0-9]*\).*/\1/p' "$log.time")
        [[ -n $signal ]] && status=$((128 + signal))
        # 124: stopped at the time limit. 137: SIGKILL -- timeout's, a minute
        # past the limit, or, before it, the kernel's out-of-memory killer.
        case $status in
            0) status=ok ;;
            124) status=timeout ;;
            137) status=$(awk -v w="$wall" -v t="$TIMEOUT" 'BEGIN { print (w >= t ? "timeout" : "oom") }') ;;
            *) status=error ;;
        esac
        # A killed run's output is partial; cooler's external sort also leaves
        # its temporary coolers behind.
        [[ $output != - && $status != ok ]] && rm -rf "$output"
        [[ $output != - ]] && rm -f "$(dirname "$output")"/tmp*.multi.cool
        bytes=0
        [[ $output != - && -e $output ]] && bytes=$(du -sb "$output" | cut -f1)
        printf '%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$section" "$dataset" "$tool" "$task" "$rep" "$status" "$wall" \
            "$(awk -v u="$user" -v s="$sys" 'BEGIN { print u + s }')" "$user" "$sys" \
            "$(awk -v k="$rss" 'BEGIN { printf "%.0f", k / 1024 }')" \
            "$(awk -v b="$bytes" 'BEGIN { printf "%.1f", b / 1e6 }')" \
            "$started" "$CONTEXT" "$command" >>"$TIMINGS"
        printf '%-44s run %d/%d  %-7s %8.2f s  %6.0f MB\n' "$key" "$rep" "$n_reps" \
            "$status" "$wall" "$(awk -v k="$rss" 'BEGIN { print k / 1024 }')"
        [[ $status == error ]] && sed 's/^/    | /' "$log.err" | tail -n 15
        [[ $status == ok ]] || break # a timeout, an OOM or a failure is not repeated
    done
    # Outputs under $SCRATCH are only timed, not used: drop them.
    [[ $output == "$SCRATCH"/* ]] && rm -rf "$output" "$output".* "${output%.*}.ps.tsv"
    return 0
}

# Record, without running it, a command that is not attempted, and why.
skip() {
    local section=$1 dataset=$2 tool=$3 task=$4 reason=$5
    read -r n_done _ <<<"$(_done "$section/$dataset/$tool/$task")"
    ((n_done > 0)) && return 0
    printf '%s\t%s\t%s\t%s\t1\tskipped\t\t\t\t\t\t\t%s\t%s\t%s\n' "$section" "$dataset" \
        "$tool" "$task" "$(date -u +%Y-%m-%dT%H:%M:%S)" "$CONTEXT" "$reason" >>"$TIMINGS"
    printf '%-44s skipped: %s\n' "$section/$dataset/$tool/$task" "$reason"
}

# Whether a command's last run in the table did not finish: it reached the
# time limit (timeout) or ran out of memory (oom).
unfinished() { [[ $(_done "$1/$2/$3/$4" | cut -d' ' -f2) =~ ^(timeout|oom)$ ]]; }

# Read files once, untimed, so that no run is the one paying the disk.
warm() { cat "$@" >/dev/null; }
