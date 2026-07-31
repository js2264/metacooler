#! /usr/bin/env bash
set -e

# --------------------------------------------------------------------------
# Hash-comparison helpers.
#
# These blocks used to just print `md5sum a b` and leave it to a human to
# notice two hashes had drifted apart — which meant a real regression could
# scroll past in CI output and the script would still exit 0. They are now
# assertions: mismatches are counted and the script exits non-zero at the end.
# Failures are accumulated rather than fatal, so one broken invariant does not
# hide the others.
# --------------------------------------------------------------------------
FAILURES=0

assert_identical() {  # assert_identical <label> <reference> <file>...
    local label="$1"; shift
    local ref="$1"; shift
    local f
    for f in "$@"; do
        if cmp -s "$ref" "$f"; then
            printf '  OK    %s: %s == %s\n' "$label" "$(basename "$f")" "$(basename "$ref")"
        else
            printf '  FAIL  %s: %s != %s\n' "$label" "$(basename "$f")" "$(basename "$ref")"
            printf '        first differing lines:\n'
            diff "$ref" "$f" | head -n 6 | sed 's/^/        /' || true
            FAILURES=$((FAILURES + 1))
        fi
    done
}

assert_differs() {  # assert_differs <label> <file_a> <file_b>
    local label="$1" a="$2" b="$3"
    if cmp -s "$a" "$b"; then
        printf '  FAIL  %s: %s == %s but they must differ\n' "$label" "$(basename "$a")" "$(basename "$b")"
        FAILURES=$((FAILURES + 1))
    else
        printf '  OK    %s: %s != %s\n' "$label" "$(basename "$a")" "$(basename "$b")"
    fi
}

mkdir -p outs
cd outs
rm -rf *

# ================================================================================
## Create metacool from pairs
# ================================================================================

for i in 1 2; do
    metacooler create pairs --contigs ../contigs_table_${i}.tsv --genomes ../mag_table_${i}.tsv --pairs ../contacts_${i}_rep1.pairs.gz --bins 1000 --output assembly_${i}_rep1_1kb.metacool
    metacooler create pairs --contigs ../contigs_table_${i}.tsv --genomes ../mag_table_${i}.tsv --pairs ../contacts_${i}_rep2.pairs.gz --bins 1000 --output assembly_${i}_rep2_1kb.metacool
    metacooler create pairs --contigs ../contigs_table_${i}.tsv --genomes ../mag_table_${i}.tsv --pairs ../contacts_${i}_rep1.pairs.gz --bins ../assembly_${i}_DpnII-digested.bed --output assembly_${i}_rep1_DpnII.metacool
    metacooler create pairs --contigs ../contigs_table_${i}.tsv --genomes ../mag_table_${i}.tsv --pairs ../contacts_${i}_rep2.pairs.gz --bins ../assembly_${i}_DpnII-digested.bed --output assembly_${i}_rep2_DpnII.metacool
done
metacooler create pairs --force --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --output assembly_1_rep1_1kb.metacool
metacooler create pairs --force --drop-unassigned --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --output assembly_1_rep1_1kb_dropped.metacool
metacooler dump -t bins assembly_1_rep1_1kb.metacool > dump_bins_before_drop.tsv
metacooler dump -t bins assembly_1_rep1_1kb_dropped.metacool > dump_bins_after_drop.tsv
git diff --no-index --compact-summary dump_bins_before_drop.tsv dump_bins_after_drop.tsv || true  ## SHOULD BE DIFFERENT, since the dropped metacooler file should have fewer bins than the original
assert_differs "--drop-unassigned removes bins" dump_bins_before_drop.tsv dump_bins_after_drop.tsv

## --contig-order fasta preserves the assembly order. It now REQUIRES a source
## for that order: --fasta, or a pairs file with #chromsize: header lines (which
## ../contacts_1_rep1.pairs.gz has). It used to fall back to alphabetical order
## silently while still recording contig-ordering="fasta".
metacooler create pairs --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --contig-order fasta --output assembly_1_rep1_1kb_same-order.metacool
cooler dump -t pixels ../assembly_1_rep1_1kb.cool > dump_pixels_cool
metacooler dump -t pixels assembly_1_rep1_1kb_same-order.metacool > dump_pixels_metacool
assert_identical "create pairs --contig-order fasta == cooler cload pairs" dump_pixels_cool dump_pixels_metacool

## -q suppresses INFO logging (warnings/errors still shown).
metacooler create pairs -q --force --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --output assembly_1_rep1_1kb_quiet.metacool

## --metadata may not carry a `commands` key: it is reserved for the provenance
## record metacooler writes itself (expected error).
printf 'commands:\n  - "not allowed"\n' > bad_metadata.yaml
metacooler create pairs --force --contigs ../contigs_table_1.tsv --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --metadata bad_metadata.yaml --output bad_metadata.metacool || echo "  (expected: 'commands' is reserved in --metadata)"

# ================================================================================
## Create metacool from cooler
# ================================================================================
metacooler create cool --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --cool ../assembly_1_rep1_1kb.cool --output assembly_1_rep1_1kb_from-cool.metacool
metacooler create cool --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --cool ../assembly_1_rep1_DpnII.cool --output assembly_1_rep1_DpnII_from-cool.metacool

cooler info ../assembly_1_rep1_1kb.cool
metacooler info assembly_1_rep1_1kb_from-cool.metacool
cooler info ../assembly_1_rep1_DpnII.cool
metacooler info assembly_1_rep1_DpnII_from-cool.metacool

# ================================================================================
## Stack / merge
# ================================================================================

metacooler stack --output assembly_1_DpnII.mmetacool assembly_1_rep1_DpnII.metacool assembly_1_rep2_DpnII.metacool
metacooler stack --output assembly_1kb.mmetacool assembly_1_rep1_1kb.metacool assembly_1_rep2_1kb.metacool assembly_2_rep1_1kb.metacool assembly_2_rep2_1kb.metacool assembly_1_rep1_DpnII.metacool assembly_1_rep2_DpnII.metacool assembly_2_rep1_DpnII.metacool assembly_2_rep2_DpnII.metacool
metacooler merge --output assembly_2_1kb.metacool assembly_2_rep1_1kb.metacool assembly_2_rep2_1kb.metacool
metacooler merge --output assembly_2_DpnII.metacool assembly_2_rep1_DpnII.metacool assembly_2_rep2_DpnII.metacool

# ================================================================================
## info/ls/tree/attrs
# ================================================================================

metacooler info assembly_1_rep1_1kb.metacool
metacooler tree assembly_1_rep1_1kb.metacool
metacooler attrs assembly_1_rep1_1kb.metacool
metacooler attrs assembly_1_DpnII.mmetacool
metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler info assembly_1kb.mmetacool::/entries/assembly_1_rep1_1kb
for file in `metacooler ls assembly_1kb.mmetacool` ; do metacooler info $file ; done

for file in `ls ../assembly_*.cool`; do metacooler ls $file; done
for file in `ls ../assembly_*.cool`; do metacooler tree $file; done
for file in `ls ../assembly_*.mcool`; do metacooler ls $file; done
for file in `ls ../assembly_*.mcool`; do metacooler tree $file; done

# ================================================================================
## all metadata subcommands
# ================================================================================

metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler metadata add --key this --value istest assembly_1_rep1_1kb.metacool
metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler metadata load ../metadata.yaml assembly_1_rep1_1kb.metacool
metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler metadata add --key this2 --value istest2 assembly_1_rep1_1kb.metacool
metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler metadata remove --key sequencing assembly_1_rep1_1kb.metacool
metacooler metadata list assembly_1_rep1_1kb.metacool
metacooler metadata get --key assembly assembly_1_rep1_1kb.metacool
metacooler metadata get --key commands assembly_1_rep1_1kb.metacool

## Every metadata subcommand needs a SINGLE map: a collection root is rejected
## (add/get/remove/load used to write to the mmetacool root, invisibly from
## every entry). Expected errors.
metacooler metadata list assembly_1kb.mmetacool || echo "  (expected: mmetacool root, specify an entry URI)"
metacooler metadata add --key x --value y assembly_1kb.mmetacool || echo "  (expected: mmetacool root, specify an entry URI)"
metacooler metadata get --key x assembly_1kb.mmetacool || echo "  (expected: mmetacool root, specify an entry URI)"
metacooler metadata load ../metadata.yaml assembly_1kb.mmetacool || echo "  (expected: mmetacool root, specify an entry URI)"
## …but an entry URI works.
metacooler metadata list assembly_1kb.mmetacool::entries/assembly_1_rep1_1kb

# ================================================================================
## ls/cp/mv/rm
# ================================================================================

metacooler ls assembly_1_DpnII.mmetacool
metacooler cp assembly_1kb.mmetacool::entries/assembly_2_rep2_1kb assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2_1kb
metacooler mv assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2_1kb assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2
metacooler rm -f assembly_1_DpnII.mmetacool::entries/assembly_1_rep1_DpnII
metacooler ls assembly_1_DpnII.mmetacool

## Answering 'y' to the overwrite prompt now really overwrites (it used to fall
## through and fail; -f was the only way).
metacooler cp assembly_1kb.mmetacool::entries/assembly_2_rep1_1kb assembly_1_DpnII.mmetacool::entries/tmp_entry
yes | metacooler cp assembly_1kb.mmetacool::entries/assembly_2_rep2_1kb assembly_1_DpnII.mmetacool::entries/tmp_entry
## Declining a confirmation prompt now exits NON-ZERO (it used to exit 0, so
## `metacooler rm … && echo ok` printed ok after you said no). Expected error.
yes n | metacooler rm assembly_1_DpnII.mmetacool::entries/tmp_entry || echo "  (expected: declining the prompt aborts with a non-zero exit)"
metacooler rm -f assembly_1_DpnII.mmetacool::entries/tmp_entry
metacooler ls assembly_1_DpnII.mmetacool

## cp/mv/rm require an ENTRY URI, not a bare collection path (expected error).
metacooler rm assembly_1_DpnII.mmetacool || echo "  (expected: entry URI required, e.g. file.mmetacool::entries/name)"

## `attrs` accepts a bare collection, but the group may only be given ONCE:
## a URI group plus a GROUP_PATH argument is rejected (expected error).
metacooler attrs assembly_1kb.mmetacool::entries/assembly_1_rep1_1kb /genomes || echo "  (expected: group given twice)"
metacooler attrs assembly_1_rep1_1kb.metacool /genomes

# ================================================================================
## cool -> metacool -> cool round trip
# ================================================================================
# these round trips should be lossless for the stored tables, ie the TABLES in the original cool file
# should be identical to the corresponding TABLES in the flattened cool file after
# going through metacool. The metacool file itself will not be identical to the
# original cool file, since it contains additional metadata and a different
# internal structure.
# Also, the final cool file after the round trip will not be exatcly identical
#  to the original cool file, because the timestamp in the header will be
# different, but the order of rows in each table HAS to be identical when using
# `--contig-order cool`

metacooler create cool --force --contig-order cool --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --output assembly_1_rep1_1kb_from_cool.metacool --cool ../assembly_1_rep1_1kb.cool
metacooler flatten --force assembly_1_rep1_1kb_from_cool.metacool -o assembly_1_rep1_1kb_from_cool_flat.cool
cooler info ../assembly_1_rep1_1kb.cool
cooler info assembly_1_rep1_1kb_from_cool_flat.cool

for table in pixels bins chroms; do
    cooler dump -t ${table} ../assembly_1_rep1_1kb.cool > dump_${table}_cool_original
    metacooler dump -t ${table} assembly_1_rep1_1kb_from_cool.metacool > dump_${table}_metacool
    cooler dump -t ${table} assembly_1_rep1_1kb_from_cool_flat.cool > dump_${table}_cool_round_trip
    git diff --no-index dump_${table}_cool_original dump_${table}_cool_round_trip
done
metacooler dump -t genomes assembly_1_rep1_1kb_from_cool.metacool > dump_genomes_metacool
metacooler dump -t contigs assembly_1_rep1_1kb_from_cool.metacool > dump_contigs_metacool
metacooler dump -t indexes assembly_1_rep1_1kb_from_cool.metacool > dump_indexes_metacool
metacooler dump -t all assembly_1_rep1_1kb_from_cool.metacool --prefix dump_all_

## --table contigs resolves the genome reference to the MAG NAME (not an integer id):
## the third column should read MAG names / __unbinned__, never 0,1,2,...
metacooler dump -t contigs -H assembly_1_rep1_1kb_from_cool.metacool | head -n 5

## --join annotates pixel dumps with contig names + coordinates (cooler --join style).
## Long: header is `contig1 start1 end1 contig2 start2 end2 count`.
metacooler dump -t pixels --join -H assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_join_metacool
head -n 5 dump_pixels_join_metacool
## For reference, cooler's own --join on the source cool (columns chrom1/2 == contig1/2).
cooler dump -t pixels --join -H ../assembly_1_rep1_1kb.cool > dump_pixels_join_cool
head -n 5 dump_pixels_join_cool

## Informational listing (the pairwise invariants are asserted above and by
## the git diff --no-index round-trip checks).
md5sum dump_*

## dump no longer clobbers -o silently: a second run needs -f/--force.
metacooler dump -t genomes -H -o dump_genomes_force.tsv assembly_1_rep1_1kb_from_cool.metacool
metacooler dump -t genomes -H -o dump_genomes_force.tsv assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: output exists, use --force)"
metacooler dump -t genomes -H -o dump_genomes_force.tsv --force -q assembly_1_rep1_1kb_from_cool.metacool

## dump rejects option/table combinations it used to accept and silently ignore.
metacooler dump -t bins --join assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: --join only applies to the pixels table)"
metacooler dump -t genomes --weight count assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: --weight only applies to the pixels table)"
metacooler dump -t indexes --genomes MAG_05 assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: --genomes cannot be used with --table indexes)"
metacooler dump -t all --prefix p_ -o out.tsv assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: -o cannot be used with --table all)"
metacooler dump -t bins --prefix p_ assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: --prefix only applies to --table all)"
## An explicitly EMPTY filter is rejected rather than silently disabled (an
## unset shell variable used to dump everything).
metacooler dump -t bins --genomes '' assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: empty --genomes value)"
## A bare .mmetacool path is rejected by every command needing a single map.
metacooler dump -t genomes assembly_1kb.mmetacool || echo "  (expected: mmetacool collection, specify an entry URI)"
metacooler info assembly_1kb.mmetacool || echo "  (expected: mmetacool collection, specify an entry URI)"

## Piping into a reader that stops early closes the pipe under dump. That used to
## print `Error: [Errno 32] Broken pipe` AFTER the rows, reading as a failure; it
## is now silent, like any other Unix tool. Nothing should appear on stderr here.
metacooler dump -t pixels assembly_1_rep1_1kb_from_cool.metacool 2> dump_head_stderr | head -n 3
if [ -s dump_head_stderr ]; then
    printf '  FAIL  dump | head wrote to stderr:\n'
    sed 's/^/        /' dump_head_stderr
    FAILURES=$((FAILURES + 1))
else
    printf '  OK    dump | head: no error reported\n'
fi

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_1_000:5000-15000 > dump_pixels_range_cool_original
metacooler dump -t pixels -r contig_1_000:5000-15000 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_1_000:5000-15000 > dump_pixels_range_cool_round_trip
assert_identical "dump -r contig:start-end" dump_pixels_range_cool_original dump_pixels_range_metacool dump_pixels_range_cool_round_trip

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_1_000:5000-15000  -r2 contig_1_001:0-10000 > dump_pixels_range2_cool_original
metacooler dump -t pixels -r contig_1_000:5000-15000  -r2 contig_1_001:0-10000 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range2_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_1_000:5000-15000  -r2 contig_1_001:0-10000 > dump_pixels_range2_cool_round_trip
assert_identical "dump -r/-r2 rectangular range" dump_pixels_range2_cool_original dump_pixels_range2_metacool dump_pixels_range2_cool_round_trip

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_1_000 > dump_pixels_range0_cool_original
metacooler dump -t pixels -r contig_1_000 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range0_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_1_000 > dump_pixels_range0_cool_round_trip
assert_identical "dump -r whole-contig range" dump_pixels_range0_cool_original dump_pixels_range0_metacool dump_pixels_range0_cool_round_trip

# ================================================================================
## extract
# ================================================================================

metacooler extract --force --output MAG_05.cool --genomes MAG_05 assembly_1_rep1_1kb.metacool
metacooler extract --force --output MAG_05-06.cool --genomes MAG_05,MAG_06 assembly_1_rep1_1kb.metacool
metacooler extract --force --output MAGs --split --genomes MAG_05,MAG_06 assembly_1_rep1_1kb.metacool
metacooler extract --force --output contig_1_000.cool --contigs contig_1_000 assembly_1_rep1_1kb.metacool
metacooler extract --force --output contig_1_000-1_001.cool --contigs contig_1_000,contig_1_001 assembly_1_rep1_1kb.metacool
## --split is genome-only: combining it with --contigs is rejected (expected error).
metacooler extract --force --output contigs.cool --split --contigs contig_1_000,contig_1_001 assembly_1_rep1_1kb.metacool || echo "  (expected: --contigs cannot be combined with --split)"

metacooler attrs assembly_1_DpnII.mmetacool

# ================================================================================
## network
# ================================================================================

metacooler network assembly_1_rep1_1kb.metacool --level contig -v -o network_contig-level.tsv
metacooler network assembly_1_rep1_1kb.metacool --level genome -v -o network_genome-level.tsv
metacooler network assembly_1_rep1_1kb.metacool --level bin -v -o network_bins-level.tsv
## network no longer clobbers -o silently either.
metacooler network assembly_1_rep1_1kb.metacool --level genome -o network_genome-level.tsv || echo "  (expected: output exists, use --force)"
metacooler network assembly_1_rep1_1kb.metacool --level genome -q --force --header --chunksize 5000 -o network_genome-level.tsv
head -n 3 network_genome-level.tsv
metacooler dump --table pixels assembly_1_rep1_1kb.metacool > dump_pixels.tsv
wc -l network_bins-level.tsv
head network_bins-level.tsv -n 20
head dump_pixels.tsv -n 20
assert_identical "network --level bin == dump -t pixels" dump_pixels.tsv network_bins-level.tsv

# ================================================================================
## aggregate
# ================================================================================

metacooler aggregate --level contig --output aggregate_contig.cool assembly_1_rep1_1kb.metacool
metacooler aggregate --level genome --output aggregate_genome.cool assembly_1_rep1_1kb.metacool
metacooler aggregate --level genome --force -q --chunksize 5000 --output aggregate_genome.cool assembly_1_rep1_1kb.metacool
cooler info aggregate_genome.cool

# ================================================================================
## binnify
# ================================================================================

## Fixed-bin source: the target MUST be an exact multiple of the source bin size.
## 5000 = 5 x 1000 -> OK.
metacooler binnify --binsize 5000 --output assembly_1kb_5kb.metacool assembly_1_rep1_1kb.metacool
metacooler info assembly_1_rep1_1kb.metacool
metacooler info assembly_1kb_5kb.metacool

## Variable-bin (DpnII fragment) source: no exact multiple exists, so the MEDIAN
## source bin width (~172 bp here) is the effective resolution; 5000 is comfortably
## coarser than that AND than the widest fragment, so no warning is emitted.
metacooler binnify --binsize 5000 --output assembly_1kb_5kb_2.metacool assembly_1_rep1_DpnII.metacool
metacooler info assembly_1_rep1_DpnII.metacool
metacooler info assembly_1kb_5kb_2.metacool

metacooler binnify --binsize 5000 --level genome --output assembly_1kb_5kb_3.metacool assembly_1_rep1_DpnII.metacool
metacooler info assembly_1kb_5kb_3.metacool

## --chunksize bounds peak memory of the pixel remapping pass; -q silences INFO.
metacooler binnify --binsize 2000 --chunksize 5000 -q --output assembly_1kb_2kb_chunked.metacool assembly_1_rep1_1kb.metacool
metacooler info assembly_1kb_2kb_chunked.metacool

## A NON-MULTIPLE target is now rejected. It used to be accepted and produced a
## file that advertised a uniform 3500 bp bin size while packing an oscillating
## 3-or-4 source bins per output bin (a built-in +/-14% checkerboard).
metacooler binnify --binsize 3500 --output assembly_1kb_3500.metacool assembly_1_rep1_1kb.metacool || echo "  (expected: 3500 is not an exact multiple of the 1000 bp source binsize)"
## A target FINER than (or equal to) the source is rejected too.
metacooler binnify --binsize 500 --output assembly_1kb_500.metacool assembly_1_rep1_1kb.metacool || echo "  (expected: new binsize must be strictly larger than the source)"
## And for a variable-bin source, a target at/below the median bin width.
metacooler binnify --binsize 50 --output assembly_DpnII_50.metacool assembly_1_rep1_DpnII.metacool || echo "  (expected: below the median source bin width of the variable-bin source)"
## Overwriting an existing output needs --force.
metacooler binnify --binsize 5000 --output assembly_1kb_5kb.metacool assembly_1_rep1_1kb.metacool || echo "  (expected: output exists, use --force)"
metacooler binnify --binsize 5000 --force --output assembly_1kb_5kb.metacool assembly_1_rep1_1kb.metacool

# ================================================================================
## zoomify
# ================================================================================

## zoomify a fixed-bin metacool into a multi-resolution .mmetacool.
## --levels N doubles the bin size per level (1000 -> 2000 -> 4000); entries are
## addressed as file.mmetacool::entries/<binsize>.
metacooler zoomify --levels 3 --output assembly_1_rep1_zoom.mmetacool assembly_1_rep1_1kb.metacool
metacooler ls assembly_1_rep1_zoom.mmetacool
metacooler info assembly_1_rep1_zoom.mmetacool::entries/1000
metacooler info assembly_1_rep1_zoom.mmetacool::entries/2000

## cooler-style explicit resolutions (--resolutions takes precedence over --levels).
## NOTE: --resolutions uses cooler's ladder-token grammar (1000N / 1000B / 4DN /
## explicit list), deliberately NOT the plain comma list taken by --genomes.
metacooler zoomify --resolutions 1000,2000,5000 --output assembly_1_rep1_zoom_r.mmetacool assembly_1_rep1_1kb.metacool
metacooler ls assembly_1_rep1_zoom_r.mmetacool

## --temp-dir puts the intermediate per-level files on a chosen scratch dir.
mkdir -p zoom_tmp
metacooler zoomify --levels 3 --temp-dir zoom_tmp -q --output assembly_1_rep1_zoom_tmp.mmetacool assembly_1_rep1_1kb.metacool
metacooler ls assembly_1_rep1_zoom_tmp.mmetacool
## …and cleans up after itself: nothing should be left behind.
ls -A zoom_tmp

## Generation STOPS EARLY instead of emitting duplicate levels: once the bin size
## exceeds the longest contig every contig collapses to one bin, so --levels 12
## yields fewer than 12 entries.
metacooler zoomify --levels 12 -q --output assembly_1_rep1_zoom_deep.mmetacool assembly_1_rep1_1kb.metacool
metacooler ls assembly_1_rep1_zoom_deep.mmetacool | wc -l

## A ladder that does not divide the base bin size is rejected: no level would
## nest in the base binning (expected error).
metacooler zoomify --resolutions 1500,3000 --output assembly_1_rep1_zoom_bad.mmetacool assembly_1_rep1_1kb.metacool || echo "  (expected: 1500 is not a multiple of the 1000 bp base binsize)"
## Variable-bin input cannot be zoomified at all (expected error).
metacooler zoomify --levels 3 --output assembly_1_rep1_zoom_dpnii.mmetacool assembly_1_rep1_DpnII.metacool || echo "  (expected: zoomify needs a fixed-bin input)"

## genome-level binnify ligates each MAG's contigs, then tiles (output contigs = MAGs).
metacooler binnify --binsize 5000 --level genome --output assembly_1_rep1_genome5kb.metacool assembly_1_rep1_1kb.metacool
metacooler dump -t contigs -H assembly_1_rep1_genome5kb.metacool | head -5

# ================================================================================
## show
# ================================================================================

## show now takes its input as a POSITIONAL path like every other command.
## -f used to mean --file; it now means --force.
metacooler show assembly_1_rep1_1kb.metacool --level genome -o show_genome.png
metacooler show assembly_1_rep1_1kb.metacool --level contig -o show_contig.png
metacooler show assembly_1_rep1_1kb.metacool --level bin --genome MAG_05 -o show_bin_MAG_05.png
## --genomes takes exactly TWO space-separated names (nargs=2), not a comma list.
metacooler show assembly_1_rep1_1kb.metacool --level bin --genomes MAG_05 MAG_06 -o show_pair.png
metacooler show assembly_1_rep1_1kb.metacool -r contig_1_000:0-10000 --transform none --cmap Reds --dpi 100 -o show_range.png
## --genome / --genomes are honoured at EVERY level (they used to be validated
## then ignored at the default --level genome).
metacooler show assembly_1_rep1_1kb.metacool --level contig --genome MAG_05 -q -o show_contig_MAG_05.png
## Re-rendering over an existing image needs -f/--force.
metacooler show assembly_1_rep1_1kb.metacool --level genome -o show_genome.png || echo "  (expected: output exists, use --force)"
metacooler show assembly_1_rep1_1kb.metacool --level genome -f -o show_genome.png
ls -1 show_*.png

# ================================================================================
## Summary
# ================================================================================
if [ "$FAILURES" -ne 0 ]; then
    printf '\nFAILED: %d hash comparison(s) did not hold.\n' "$FAILURES"
    exit 1
fi
printf '\nAll hash comparisons held.\n'
