# shellcheck shell=bash
# Sourced by every bash cell of benchmarks.qmd: settings, then `bench`.
#
#   bench SECTION DATASET TOOL TASK OUTPUT -- COMMAND [ARG...]
#   EVICT="FILE..." bench SECTION DATASET TOOL TASK OUTPUT -- COMMAND [ARG...]
#   both "FILE..." SECTION DATASET TOOL TASK OUTPUT -- COMMAND [ARG...]
#
# `bench` runs COMMAND as a fresh process, each run under GNU time and
# `timeout`, after deleting OUTPUT (the file it writes, '-' for none) and its
# sidecars; and appends one row per run to $RESULTS/timings.tsv, with its
# status: ok, timeout (stopped at $TIMEOUT s), oom (killed by the kernel for
# want of memory), error, or skipped.
#
# Cache. A run is warm by default: what the command reads may be in the page
# cache, as it is right after it was written or read. With EVICT, it is cold:
# the FILEs (what the command reads) are evicted from the page cache before
# every run (pagecache.py), so their bytes come from the disk. `both` times a
# command warm (its FILEs read just before), then cold.
#
# Runs. At least $MIN_REPS, then more, up to $MAX_REPS, while the runs so far
# add up to under $BUDGET s. A command whose run takes longer than $LONG s
# runs once. BENCH_REPS=N runs every command N times instead. A timeout, an
# OOM or a failure is not repeated.
#
# Peak memory is GNU time's: the largest resident set among the command's
# processes; for a pipeline whose stages are each prefixed with $STAGE, the
# sum of the stages' own peaks, since they run side by side. A pipeline that
# does not finish has no peak (left empty): its shell is killed before it
# reaps the stages, so their memory never reaches GNU time.
#
# A command already in the table, with as many runs as it would get and its
# output still there, is not run again (BENCH_REDO=1 forces it), so a page can
# be re-rendered without re-timing.

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
MIN_REPS=${BENCH_MIN_REPS:-3}          # runs per command, at least
MAX_REPS=${BENCH_MAX_REPS:-5}          # ... and at most
BUDGET=${BENCH_BUDGET:-120}            # seconds of runs past which no more than MIN_REPS
LONG=${BENCH_LONG:-600}                # seconds: a command whose run takes longer runs once
COOLER_VERSION=0.12.1                  # the cooler every timing is against: uv.lock's
ZOOMS=${BENCH_RESOLUTIONS:-$BINSIZE,$((2 * BINSIZE)),$((5 * BINSIZE)),$((10 * BINSIZE))}

# Prefix for each stage of a timed pipeline: its peak memory is added up with
# the others'. Single-quoted: the pipeline's own shell expands $BENCH_STAGES.
STAGE='/usr/bin/time -f %M -a -o "$BENCH_STAGES"'

mkdir -p "$DATA" "$MAPS" "$QUERIES" "$SCRATCH" "$LOGS" "$RESULTS"

COLUMNS_TSV="section	dataset	tool	task	cache	rep	status	wall_s	cpu_s	user_s	sys_s	peak_mb	output_mb	started	metacooler	cooler	cpu	cores	memory_gb	command"
# A table an older bench.sh wrote has other columns: keep it aside, start anew.
if [[ -s $TIMINGS && $(head -n 1 "$TIMINGS") != "$COLUMNS_TSV" ]]; then
    aside=$WORK/timings.$(date -u +%Y%m%dT%H%M%S).tsv
    mv "$TIMINGS" "$aside"
    echo "timings.tsv has the columns of an older bench.sh: moved to $aside; timing anew" >&2
fi
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

# A command's rows so far: how many, the seconds its finished runs took in
# all, the seconds its last run took, and that run's status ("none" yet).
_done() {
    awk -F'\t' -v k="$1" '$1"/"$2"/"$3"/"$4"/"$5 == k { n++; s = $7; l = $8; if ($7 == "ok") t += $8 }
        END { print n + 0, t + 0, l + 0, (n ? s : "none") }' "$TIMINGS"
}

# Whether a command that ran N times, for TOTAL s, the last for LAST s, runs again.
_more() {
    if [[ -n ${BENCH_REPS:-} ]]; then
        (($1 < BENCH_REPS))
        return
    fi
    awk -v n="$1" -v t="$2" -v l="$3" -v lo="$MIN_REPS" -v hi="$MAX_REPS" -v b="$BUDGET" -v long="$LONG" \
        'BEGIN { exit !(n == 0 || (n < hi && l <= long && (n < lo || t < b))) }'
}

# Record, without running it, a command that is not attempted, and why.
_skip() {
    local section=$1 dataset=$2 tool=$3 task=$4 cache=$5 reason=$6 n
    read -r n _ <<<"$(_done "$section/$dataset/$tool/$task/$cache")"
    ((n > 0)) && return 0
    printf '%s\t%s\t%s\t%s\t%s\t1\tskipped\t\t\t\t\t\t\t%s\t%s\t%s\n' "$section" "$dataset" "$tool" \
        "$task" "$cache" "$(date -u +%Y-%m-%dT%H:%M:%S)" "$CONTEXT" "$reason" >>"$TIMINGS"
    printf '%-52s skipped: %s\n' "$section/$dataset/$tool/$task/$cache" "$reason"
}
skip() { _skip "$1" "$2" "$3" "$4" warm "$5"; }

bench() {
    local section=$1 dataset=$2 tool=$3 task=$4 output=$5
    shift 5
    [[ ${1-} == -- ]] && shift
    local cache=warm evict=()
    if [[ -n ${EVICT:-} ]]; then
        cache=cold
        read -ra evict <<<"$EVICT"
    fi
    local key="$section/$dataset/$tool/$task/$cache" n total last_wall last remake=0
    if [[ ${BENCH_REDO:-0} == 1 ]]; then
        awk -F'\t' -v k="$key" '$1"/"$2"/"$3"/"$4"/"$5 != k' "$TIMINGS" >"$TIMINGS.tmp"
        mv "$TIMINGS.tmp" "$TIMINGS"
    fi
    read -r n total last_wall last <<<"$(_done "$key")"
    if [[ $last =~ ^(timeout|oom|skipped)$ ]]; then # final: it would only fail again
        printf '%-52s did not finish before (%s): not run again\n' "$key" "$last"
        return 0
    fi
    if ! _more "$n" "$total" "$last_wall"; then
        if [[ $output == - || $output == "$SCRATCH"/* || -e $output || $last != ok ]]; then
            printf '%-52s already timed (%d runs, last %s)\n' "$key" "$n" "$last"
            return 0
        fi
        remake=1 # its output is gone: one more run to make it again
    fi
    local log="$LOGS/${key//\//_}" command="$*" rep status signal wall user sys rss stages peak bytes started
    local staged=0
    [[ $* == *"$STAGE"* ]] && staged=1
    command=${command//"$STAGE "/} # the wrapper measures the command; it is not part of it
    command=${command//$'\t'/ }
    command=${command//$'\n'/ }
    while ((remake)) || _more "$n" "$total" "$last_wall"; do
        remake=0
        rep=$((n + 1))
        if ((${#evict[@]})) && ! python "$HERE/pagecache.py" evict "${evict[@]}" 2>"$log.evict"; then
            _skip "$section" "$dataset" "$tool" "$task" "$cache" \
                "cannot evict its inputs from the page cache: $(head -n 1 "$log.evict")"
            return 0
        fi
        # The output, and the P(s) table `create pairs` writes beside it.
        [[ $output != - ]] && rm -rf "$output" "$output".* "${output%.*}.ps.tsv"
        : >"$log.stages"
        started=$(date -u +%Y-%m-%dT%H:%M:%S)
        BENCH_STAGES=$log.stages /usr/bin/time -f '%e %U %S %M %x' -o "$log.time" \
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
        # A pipeline's stages run side by side: its peak is theirs added up.
        stages=$(awk '/^[0-9]+$/ { s += $1 } END { print s + 0 }' "$log.stages")
        ((stages > rss)) && rss=$stages
        peak=$(awk -v k="$rss" 'BEGIN { printf "%.0f", k / 1024 }')
        [[ $staged == 1 && $status != ok ]] && peak="" # unknown: see above
        # A killed run's output is partial; cooler's external sort also leaves
        # its temporary coolers behind.
        [[ $output != - && $status != ok ]] && rm -rf "$output"
        [[ $output != - ]] && rm -f "$(dirname "$output")"/tmp*.multi.cool
        bytes=0
        [[ $output != - && -e $output ]] && bytes=$(du -sb "$output" | cut -f1)
        printf '%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$section" "$dataset" "$tool" "$task" "$cache" "$rep" "$status" "$wall" \
            "$(awk -v u="$user" -v s="$sys" 'BEGIN { print u + s }')" "$user" "$sys" "$peak" \
            "$(awk -v b="$bytes" 'BEGIN { printf "%.1f", b / 1e6 }')" \
            "$started" "$CONTEXT" "$command" >>"$TIMINGS"
        printf '%-52s run %d  %-7s %8.2f s  %6s MB\n' "$key" "$rep" "$status" "$wall" "${peak:-?}"
        [[ $status == error ]] && sed 's/^/    | /' "$log.err" | tail -n 15
        [[ $status == ok ]] || break # a timeout, an OOM or a failure is not repeated
        n=$rep
        total=$(awk -v t="$total" -v w="$wall" 'BEGIN { print t + w }')
        last_wall=$wall
    done
    # Outputs under $SCRATCH are only timed, not used: drop them.
    [[ $output == "$SCRATCH"/* ]] && rm -rf "$output" "$output".* "${output%.*}.ps.tsv"
    return 0
}

# Time a command warm, its FILEs read just before, then cold, its FILEs
# evicted from the page cache before every run.
both() {
    local files=$1
    shift
    # shellcheck disable=SC2086 # FILEs are a list
    warm $files
    bench "$@"
    EVICT=$files bench "$@"
}

# Whether a command's last (warm) run in the table did not finish: it reached
# the time limit (timeout) or ran out of memory (oom).
unfinished() { [[ $(_done "$1/$2/$3/$4/warm" | cut -d' ' -f4) =~ ^(timeout|oom)$ ]]; }

# Read files once, untimed, so that no run is the one paying the disk.
warm() { cat "$@" >/dev/null; }
