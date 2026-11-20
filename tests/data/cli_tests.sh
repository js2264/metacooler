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
## Collect single-genome cools into one metacool (one genome per input)
# ================================================================================

printf 'genome_name\tsource\nA1\tassembly_1\nA2\tassembly_2\n' > collected_species.tsv
metacooler collect ../assembly_1_rep1_1kb.cool ../assembly_2_rep1_1kb.cool --names A1,A2 --genomes collected_species.tsv --output collected_1kb.metacool
metacooler info collected_1kb.metacool
metacooler dump -t genomes -H collected_1kb.metacool
metacooler dump -t contigs -H collected_1kb.metacool | head -n 5

## Under the default --contig-order cool, extracting a genome gives back its
## input: the same bins and pixels (only the contig names carry the A2__ prefix).
metacooler extract collected_1kb.metacool --genomes A2 --output collected_A2.cool
cooler dump -t pixels ../assembly_2_rep1_1kb.cool > dump_pixels_collect_original
cooler dump -t pixels collected_A2.cool > dump_pixels_collect_round_trip
assert_identical "collect -> extract round trip" dump_pixels_collect_original dump_pixels_collect_round_trip

## Reordering contigs within a genome re-sorts only what moved; the map is the same.
metacooler collect ../assembly_1_rep1_1kb.cool ../assembly_2_rep1_1kb.cool --names A1,A2 --contig-order genome-then-length --sort-strategy external --output collected_1kb_by_length.metacool
metacooler check --deep collected_1kb_by_length.metacool

## One resolution from each .mcool.
metacooler collect ../assembly_1_rep1_1kb.mcool::resolutions/5000 ../assembly_2_rep1_1kb.mcool::resolutions/5000 --output collected_5kb.metacool
metacooler info collected_5kb.metacool

## Refused inputs (expected errors).
metacooler collect ../assembly_1_rep1_1kb.cool ../assembly_2_rep1_DpnII.cool --output collected_mixed.metacool || echo "  (expected: inputs are not binned alike)"
metacooler collect ../assembly_1_rep1_1kb.mcool ../assembly_2_rep1_1kb.cool --output collected_mcool.metacool || echo "  (expected: pick one resolution of the .mcool)"
metacooler collect assembly_1_rep1_1kb.metacool ../assembly_2_rep1_1kb.cool --output collected_metacool.metacool || echo "  (expected: already a metacool file)"

# ================================================================================
## Stack / merge
# ================================================================================

metacooler stack --output assembly_1_DpnII.mmetacool assembly_1_rep1_DpnII.metacool assembly_1_rep2_DpnII.metacool
metacooler stack --output assembly_1kb.mmetacool assembly_1_rep1_1kb.metacool assembly_1_rep2_1kb.metacool assembly_2_rep1_1kb.metacool assembly_2_rep2_1kb.metacool assembly_1_rep1_DpnII.metacool assembly_1_rep2_DpnII.metacool assembly_2_rep1_DpnII.metacool assembly_2_rep2_DpnII.metacool
metacooler merge --output assembly_2_1kb.metacool assembly_2_rep1_1kb.metacool assembly_2_rep2_1kb.metacool
metacooler merge --output assembly_2_DpnII.metacool assembly_2_rep1_DpnII.metacool assembly_2_rep2_DpnII.metacool

## Both replicates' pairs files in one create pairs pass: the contacts and the
## QC block (with its qc_* columns) merge sums, plus the P(s) table merge cannot give.
metacooler create pairs -q --contigs ../contigs_table_2.tsv --genomes ../mag_table_2.tsv --pairs ../contacts_2_rep1.pairs.gz --pairs ../contacts_2_rep2.pairs.gz --bins 1000 --output assembly_2_1kb_one-pass.metacool
for table in pixels genomes contigs; do
    metacooler dump -t ${table} assembly_2_1kb.metacool > dump_${table}_merged_reps
    metacooler dump -t ${table} assembly_2_1kb_one-pass.metacool > dump_${table}_one_pass_reps
    assert_identical "create pairs --pairs rep1 --pairs rep2 == merge rep1 rep2 (${table})" dump_${table}_merged_reps dump_${table}_one_pass_reps
done
metacooler metadata get --key qc assembly_2_1kb.metacool > qc_merged_reps
metacooler metadata get --key qc assembly_2_1kb_one-pass.metacool > qc_one_pass_reps
assert_identical "create pairs --pairs rep1 --pairs rep2 == merge rep1 rep2 (QC block)" qc_merged_reps qc_one_pass_reps

# ================================================================================
## reassign
# ================================================================================

## A reassignment is a relabelling: the pixel table is identical, whichever path
## is taken, so a `dump -t pixels` of the source and of every output must match.
metacooler dump -t pixels assembly_1_rep1_1kb.metacool > reassign_pixels_source.tsv
metacooler dump -t bins assembly_1_rep1_1kb.metacool > reassign_bins_source.tsv
metacooler dump -t contigs assembly_1_rep1_1kb.metacool > reassign_contigs_source.tsv
metacooler dump -t genomes assembly_1_rep1_1kb.metacool > reassign_genomes_source.tsv

## Incremental (--missing keep, the default): only the named contigs move.
printf 'contig_name\tgenome_name\ncontig_1_000\tMAG_02\n' > reassign_one.tsv
metacooler reassign assembly_1_rep1_1kb.metacool --contigs reassign_one.tsv -o reassign_one.metacool
metacooler check reassign_one.metacool
metacooler dump -t genomes reassign_one.metacool > reassign_one_genomes.tsv
assert_differs "reassign moved a contig" reassign_genomes_source.tsv reassign_one_genomes.tsv

## --layout keep never re-orders, so the bins are untouched and the file records
## contig-ordering="custom" (it can no longer claim to be genome-grouped).
metacooler reassign assembly_1_rep1_1kb.metacool --contigs reassign_one.tsv --layout keep -o reassign_keep.metacool
metacooler check reassign_keep.metacool
metacooler attrs reassign_keep.metacool | grep contig-ordering
metacooler dump -t bins reassign_keep.metacool > reassign_keep_bins.tsv
assert_identical "reassign --layout keep leaves /bins alone" reassign_bins_source.tsv reassign_keep_bins.tsv

## Whatever the layout, no contact is created, destroyed or moved between bins.
metacooler dump -t pixels reassign_keep.metacool > reassign_keep_pixels.tsv
assert_identical "reassign preserves the pixel table" reassign_pixels_source.tsv reassign_keep_pixels.tsv

## The creation table, re-applied wholesale, is a no-op on the assignment — a
## full round trip through the re-ordering path.
metacooler reassign assembly_1_rep1_1kb.metacool --contigs ../contigs_table_1.tsv --missing unbinned -o reassign_roundtrip.metacool
metacooler dump -t contigs reassign_roundtrip.metacool > reassign_roundtrip_contigs.tsv
assert_identical "re-applying the creation table changes nothing" reassign_contigs_source.tsv reassign_roundtrip_contigs.tsv

## In place, on a copy: allowed only because this reassignment moves no contig.
cp assembly_1_rep1_1kb.metacool reassign_inplace.metacool
metacooler reassign reassign_inplace.metacool --contigs reassign_one.tsv --layout keep --inplace
metacooler check reassign_inplace.metacool

## Expected errors: a table for another assembly, and --inplace when the layout
## has to be restored.
printf 'contig_name\tgenome_name\nnot_a_contig\tMAG_02\n' > reassign_wrong.tsv
metacooler reassign assembly_1_rep1_1kb.metacool --contigs reassign_wrong.tsv -o reassign_wrong.metacool || echo "  (expected: no contig name matches this file)"
cp assembly_1_rep1_1kb.metacool reassign_inplace_bad.metacool
metacooler reassign reassign_inplace_bad.metacool --contigs reassign_one.tsv --layout genome-then-length --inplace || echo "  (expected: --inplace cannot re-order contigs)"

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
## the second column should read MAG names / __unbinned__, never 0,1,2,...
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
metacooler dump -t indexes --filter 'contigs.length > 1' assembly_1_rep1_1kb_from_cool.metacool || echo "  (expected: --filter cannot be used with --table indexes)"
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

# ================================================================================
## --filter
# ================================================================================

## A --filter at CREATE time and the same --filter at READ time must select the
## same contigs: the create path evaluates the expression against tables built
## by the very same builders the output file uses, so the two can only disagree
## if that stops being true.
metacooler create pairs -q --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv \
    --pairs ../contacts_1_rep1.pairs.gz --bins 1000 --filter 'contigs.length >= 30kb' \
    --output assembly_1_rep1_1kb_filtered.metacool
metacooler dump -t contigs -H assembly_1_rep1_1kb_filtered.metacool | cut -f1,2 | sort > filter_contigs_at_create.tsv
metacooler dump -t contigs -H --filter 'contigs.length >= 30kb' assembly_1_rep1_1kb.metacool | cut -f1,2 | sort > filter_contigs_at_read.tsv
assert_identical "--filter at create == --filter at dump" filter_contigs_at_read.tsv filter_contigs_at_create.tsv

## A filter is a genuine restriction, not a no-op.
metacooler dump -t contigs -H assembly_1_rep1_1kb.metacool | cut -f1,2 | sort > filter_contigs_unfiltered.tsv
assert_differs "--filter restricts the contigs table" filter_contigs_unfiltered.tsv filter_contigs_at_read.tsv

## --genomes and --filter INTERSECT: the result must be a proper, non-empty
## subset of what --genomes alone gives, and must never pull in a contig from
## another genome (which a union would). MAG_01's contigs are 25 kb, 22 kb and
## 13 kb, so a 20 kb threshold keeps two of the three while the file as a whole
## has longer contigs elsewhere.
metacooler dump -t contigs --genomes MAG_01 assembly_1_rep1_1kb.metacool | cut -f1 | sort > filter_mag01_all.txt
metacooler dump -t contigs --genomes MAG_01 --filter 'contigs.length >= 20kb' assembly_1_rep1_1kb.metacool | cut -f1 | sort > filter_mag01_long.txt
foreign="$(comm -13 filter_mag01_all.txt filter_mag01_long.txt)"
n_all="$(wc -l < filter_mag01_all.txt)"
n_long="$(wc -l < filter_mag01_long.txt)"
if [ -n "$foreign" ]; then
    printf '  FAIL  --genomes + --filter returned contigs outside MAG_01:\n'
    printf '%s\n' "$foreign" | sed 's/^/        /'
    FAILURES=$((FAILURES + 1))
elif [ "$n_long" -eq 0 ] || [ "$n_long" -ge "$n_all" ]; then
    printf '  FAIL  --genomes + --filter kept %s of %s contigs — not a proper subset, so the intersection is untested\n' "$n_long" "$n_all"
    FAILURES=$((FAILURES + 1))
else
    printf '  OK    --genomes + --filter intersect: %s of MAG_01'"'"'s %s contigs, none foreign\n' "$n_long" "$n_all"
fi

## extract accepts --filter as a complete selection, with no name given.
metacooler extract assembly_1_rep1_1kb.metacool --filter 'genomes.length > 100kb' -o filter_extract.cool
metacooler show assembly_1_rep1_1kb.metacool --filter 'genomes.length > 100kb' -o filter_show.png

## Malformed / impossible filters are errors, not empty output (expected errors).
metacooler dump -t genomes --filter 'length > 1Mb' assembly_1_rep1_1kb.metacool || echo "  (expected: --filter needs a genomes./contigs. table prefix)"
metacooler dump -t genomes --filter 'genomes.lenght > 1Mb' assembly_1_rep1_1kb.metacool || echo "  (expected: unknown column, with a suggestion)"
metacooler dump -t genomes --filter 'genomes.length > 1Gb' assembly_1_rep1_1kb.metacool || echo "  (expected: nothing passed --filter)"

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
## There is no singular --genome any more: a one-name --genomes IS that.
metacooler show assembly_1_rep1_1kb.metacool --level bin --genomes MAG_05 -o show_bin_MAG_05.png
metacooler show assembly_1_rep1_1kb.metacool --level genome --genome MAG_05 -o show_singular.png || echo "  (expected: no such option --genome)"
## --genomes is a comma-separated / repeatable list, and it is SYMMETRIC for any
## number of names: every name goes on BOTH axes, so the image holds all the
## intra- and inter-genome blocks between them, in the order given.
metacooler show assembly_1_rep1_1kb.metacool --level bin --genomes MAG_05,MAG_06 -o show_pair.png
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_05 --genomes MAG_06 -o show_pair_contig.png
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_05,MAG_06,MAG_07,MAG_08 -o show_multi.png
metacooler show assembly_1_rep1_1kb.metacool --level bin --genomes MAG_05,MAG_06,MAG_07 --dpi 100 -o show_multi_bin.png
metacooler show assembly_1_rep1_1kb.metacool --level genome --genomes MAG_05,MAG_06,MAG_07 -o show_multi_genome.png
## The block order follows the order the names were given, so these two differ.
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_05,MAG_06,MAG_07 --dpi 100 -o show_multi_order_a.png
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_07,MAG_05,MAG_06 --dpi 100 -o show_multi_order_b.png
assert_differs "show --genomes block order" show_multi_order_a.png show_multi_order_b.png
## --contigs is the same selector in contig space -- the one thing -r cannot do
## is put SEVERAL contigs on the same axis.
metacooler show assembly_1_rep1_1kb.metacool --level bin --contigs contig_1_000,contig_1_010 --dpi 100 -o show_contigs_bin.png
metacooler show assembly_1_rep1_1kb.metacool --level contig --contigs contig_1_000 --contigs contig_1_010,contig_1_023 -o show_contigs_contig.png
## A contig has no row of its own at --level genome, so the level is upgraded.
metacooler show assembly_1_rep1_1kb.metacool --level genome --contigs contig_1_000,contig_1_010 -v -o show_contigs_upgraded.png
## ...but the two selectors may not be combined: both put their names on BOTH
## axes, so together they have no single meaning (expected error).
metacooler show assembly_1_rep1_1kb.metacool --genomes MAG_05 --contigs contig_1_000 -o show_both.png || echo "  (expected: --genomes and --contigs are mutually exclusive)"
## A name may only appear once, and an empty value is rejected (expected errors).
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_05,MAG_06,MAG_05 -o show_dup.png || echo "  (expected: MAG_05 listed more than once)"
metacooler show assembly_1_rep1_1kb.metacool --level contig --contigs '' -o show_empty_contigs.png || echo "  (expected: empty --contigs value)"
## The old space-separated form is now a usage error: the second name would fall
## through to the positional METACOOL_PATH slot (expected error).
metacooler show assembly_1_rep1_1kb.metacool --level bin --genomes MAG_05 MAG_06 -o show_spaced.png || echo "  (expected: --genomes is a comma list, not nargs=2)"
## -r/-r2 is the ONLY asymmetric form (a different region per axis), and it ticks
## the axes in BASE PAIRS, not bin indices (a bin index starts at 0 even when the
## range does not).
metacooler show assembly_1_rep1_1kb.metacool -r contig_1_000:0-10000 --transform none --cmap Reds --dpi 100 -o show_range.png
metacooler show assembly_1_rep1_1kb.metacool -r contig_1_010:3000-9000 --dpi 100 -o show_range_offset.png
metacooler show assembly_1_rep1_1kb.metacool -r MAG_05 -r2 MAG_06 --dpi 100 -o show_range_rect.png
## --genomes / --contigs are honoured at EVERY level (they used to be validated
## then ignored at the default --level genome).
metacooler show assembly_1_rep1_1kb.metacool --level contig --genomes MAG_05 -q -o show_contig_MAG_05.png
## Re-rendering over an existing image needs -f/--force.
metacooler show assembly_1_rep1_1kb.metacool --level genome -o show_genome.png || echo "  (expected: output exists, use --force)"
metacooler show assembly_1_rep1_1kb.metacool --level genome -f -o show_genome.png
ls -1 show_*.png

# ================================================================================
## qc
# ================================================================================

## `create` already stored a block measured on the raw read positions.
metacooler metadata get assembly_1_rep1_1kb.metacool -k qc
## The per-genome columns live in /genomes and dump like any other column.
metacooler dump --table genomes --header assembly_1_rep1_1kb.metacool | head -n 3
## Recomputing from the pixel table gives the same totals with a coarser
## intra-contig split, and says so via "source": "pixels".
metacooler qc assembly_1_rep1_1kb.metacool -q
metacooler qc assembly_1_rep1_1kb.metacool --table --threshold 5kb -q
metacooler qc assembly_1_rep1_1kb.metacool --per-genome qc_per_genome.tsv -q > qc_pixels.json
head -n 3 qc_per_genome.tsv
## Every category is a share of the same total, so the four sum back to it.
python -c "
import json
b = json.load(open('qc_pixels.json'))
parts = ['intra_contig_short', 'intra_contig_long', 'inter_contig_intra_genome', 'inter_genome']
assert sum(b[k] for k in parts) == b['n_contacts'], b
print('  OK    qc: categories sum to n_contacts')
"
## -o refuses to clobber, like every other writer (expected error).
metacooler qc assembly_1_rep1_1kb.metacool -o qc_pixels.json -q || echo "  (expected: output exists, use --force)"
## --store writes the block back; the file must stay a valid metacool AND a
## valid cooler afterwards.
cp assembly_1_rep1_1kb.metacool qc_stored.metacool
metacooler qc qc_stored.metacool --store -q > /dev/null
metacooler check qc_stored.metacool
cooler info qc_stored.metacool > /dev/null && echo "  OK    qc --store: still opens as a cooler"

# ================================================================================
## ps
# ================================================================================

## `create pairs` already wrote the fine-grained, per-MAG table beside its
## output — measured on read positions, so it resolves below the 1 kb grid.
head -n 3 assembly_1_rep1_1kb.ps.tsv
## Recomputing it from the same pairs must reproduce that file byte for byte.
metacooler ps assembly_1_rep1_1kb.metacool --pairs ../contacts_1_rep1.pairs.gz \
    --by genome -o ps_recomputed.tsv -q
assert_identical "ps --pairs == create pairs sidecar" \
    assembly_1_rep1_1kb.ps.tsv ps_recomputed.tsv
## …and the pairs-derived curve reaches distances the binned one cannot.
python -c "
import pandas as pd
pairs = pd.read_csv('ps_recomputed.tsv', sep='\t')
assert pairs['dist'].min() < 1000, pairs['dist'].min()
print('  OK    ps --pairs: resolves below the bin size')
"

metacooler ps assembly_1_rep1_1kb.metacool -q | head -n 5
metacooler ps assembly_1_rep1_1kb.metacool -q -o ps_pooled.tsv --plot ps_pooled.png
metacooler ps assembly_1_rep1_1kb.metacool -q --by genome -o ps_by_genome.tsv --plot ps_by_genome.png
metacooler ps assembly_1_rep1_1kb.metacool -q --genomes MAG_05,MAG_06 --bins-per-decade 5 -o ps_two_mags.tsv
metacooler ps assembly_1_rep1_1kb.metacool -q --filter "genomes.length > 100kb" --value avg_contacts --plot ps_big.png -o ps_big.tsv
## Summing the per-genome curves returns the pooled one: same contacts, same
## possible bin pairs, just grouped differently.
python -c "
import pandas as pd
pooled = pd.read_csv('ps_pooled.tsv', sep='\t')
by_genome = pd.read_csv('ps_by_genome.tsv', sep='\t')
assert pooled['n_contacts'].sum() == by_genome['n_contacts'].sum()
assert pooled['n_bin_pairs'].sum() == by_genome['n_bin_pairs'].sum()
print('  OK    ps: --by genome partitions the pooled curve')
"
## Both outputs refuse to clobber (expected error).
metacooler ps assembly_1_rep1_1kb.metacool -q -o ps_pooled.tsv || echo "  (expected: output exists, use --force)"
ls -1 ps_*.png ps_*.tsv

# ================================================================================
## check
# ================================================================================

## The integrity report over the existing validators: one line per check instead
## of raising on the first problem.
metacooler check assembly_1_rep1_1kb.metacool
## --deep additionally scans the WHOLE pixel table and reconciles attrs['sum'].
metacooler check --deep assembly_1_rep1_1kb.metacool
## A bare collection path is accepted: the per-map checks repeat for each entry.
metacooler check assembly_1kb.mmetacool
metacooler check assembly_1kb.mmetacool::entries/assembly_1_rep1_1kb
## A plain .cool is reported as such rather than as a broken metacool, and exits
## non-zero (expected error).
metacooler check ../assembly_1_rep1_1kb.cool || echo "  (expected: not a metacool file)"
## Corrupting a declared count is caught by the structure check (expected error).
cp assembly_1_rep1_1kb.metacool check_broken.metacool
python -c "
import h5py
with h5py.File('check_broken.metacool', 'r+') as f:
    f.attrs['nnz'] = int(f.attrs['nnz']) + 1
"
metacooler check check_broken.metacool || echo "  (expected: nnz disagrees with /pixels)"

# ================================================================================
## Summary
# ================================================================================
if [ "$FAILURES" -ne 0 ]; then
    printf '\nFAILED: %d hash comparison(s) did not hold.\n' "$FAILURES"
    exit 1
fi
printf '\nAll hash comparisons held.\n'
