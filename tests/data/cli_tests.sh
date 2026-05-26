#! /usr/bin/env bash
set -e

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
git diff --no-index dump_bins_before_drop.tsv dump_bins_after_drop.tsv
md5sum dump_bins_before_drop.tsv dump_bins_after_drop.tsv

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

# ================================================================================
## ls/cp/mv/rm
# ================================================================================

metacooler ls assembly_1_DpnII.mmetacool
metacooler cp assembly_1kb.mmetacool::entries/assembly_2_rep2_1kb assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2_1kb
metacooler mv assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2_1kb assembly_1_DpnII.mmetacool::/entries/assembly_2_rep2
metacooler rm -f assembly_1_DpnII.mmetacool::entries/assembly_1_rep1_DpnII
metacooler ls assembly_1_DpnII.mmetacool

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

metacooler create cool --force --contigs ../contigs_table_1.tsv --genomes ../mag_table_1.tsv --output assembly_1_rep1_1kb_from_cool.metacool --cool ../assembly_1_rep1_1kb.cool
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

md5sum dump_*

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_0_0:5000-15000 > dump_pixels_range_cool_original
metacooler dump -t pixels -r contig_0_0:5000-15000 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_0_0:5000-15000 > dump_pixels_range_cool_round_trip
md5sum dump_pixels_range_*

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_0_0:5000-15000  -r2 contig_0_1:0-10000 > dump_pixels_range2_cool_original
metacooler dump -t pixels -r contig_0_0:5000-15000  -r2 contig_0_1:0-10000 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range2_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_0_0:5000-15000  -r2 contig_0_1:0-10000 > dump_pixels_range2_cool_round_trip
md5sum dump_pixels_range2_*

cooler dump -t pixels ../assembly_1_rep1_1kb.cool -r contig_0_0 > dump_pixels_range0_cool_original
metacooler dump -t pixels -r contig_0_0 assembly_1_rep1_1kb_from_cool.metacool > dump_pixels_range0_metacool
cooler dump -t pixels assembly_1_rep1_1kb_from_cool_flat.cool -r contig_0_0 > dump_pixels_range0_cool_round_trip
md5sum dump_pixels_range0_*

# ================================================================================
## extract
# ================================================================================

metacooler extract --force --output MAG_05.cool --genomes MAG_05 assembly_1_rep1_1kb.metacool
metacooler extract --force --output MAG_05-06.cool --genomes MAG_05,MAG_06 assembly_1_rep1_1kb.metacool
metacooler extract --force --output MAGs --split --genomes MAG_05,MAG_06 assembly_1_rep1_1kb.metacool
metacooler extract --force --output contig_0_0.cool --contigs contig_0_0 assembly_1_rep1_1kb.metacool
metacooler extract --force --output contig_0_0-1_0.cool --contigs contig_0_0,contig_1_0 assembly_1_rep1_1kb.metacool
metacooler extract --force --output contigs.cool --split --contigs contig_0_0,contig_1_0 assembly_1_rep1_1kb.metacool

metacooler attrs assembly_1_DpnII.mmetacool

# ================================================================================
## network
# ================================================================================

metacooler network assembly_1_rep1_1kb.metacool --level contig -v -o network_contig-level.tsv
metacooler network assembly_1_rep1_1kb.metacool --level genome -v -o network_genome-level.tsv
metacooler network assembly_1_rep1_1kb.metacool --level bin -v -o network_bins-level.tsv
metacooler dump --table pixels assembly_1_rep1_1kb.metacool > dump_pixels.tsv
wc -l network_bins-level.tsv
head network_bins-level.tsv -n 20
head dump_pixels.tsv -n 20
md5sum network_bins-level.tsv dump_pixels.tsv

# ================================================================================
# ================================================================================
# ================================================================================
# ================================================================================
# ================================================================================
## TO DO
# ================================================================================
# ================================================================================
# ================================================================================
# ================================================================================
# ================================================================================

## │ aggregate         Aggregate contacts to contig or genome level.                                                       │
## │ binnify            Re-bin a metacool file to a new bin size.                                                          │
## │ zoomify            Create a multi-resolution mmetacool file.                                                          │
## │ show              Visualize contact maps.                                                                             │

# ================================================================================
## binnify
# ================================================================================

metacooler binnify --binsize 5000 --output assembly_1kb_5kb.metacool assembly_1_rep1_1kb.metacool
metacooler info assembly_1_rep1_1kb.metacool
metacooler info assembly_1kb_5kb.metacool

metacooler binnify --binsize 5000 --output assembly_1kb_5kb_2.metacool assembly_1_rep1_DpnII.metacool
metacooler info assembly_1_rep1_DpnII.metacool
metacooler info assembly_1kb_5kb_2.metacool

metacooler binnify --binsize 5000 --level genomes --output assembly_1kb_5kb_3.metacool assembly_1_rep1_DpnII.metacool
metacooler info assembly_1kb_5kb_3.metacool

# ================================================================================
## zoomify
# ================================================================================

metacooler zoomify --binsize 1000,2000,5000 --output assembly_1kb_5kb.metacool assembly_1_rep1_1kb.metacool
metacooler info assembly_1_rep1_1kb.mmetacool
metacooler info assembly_1kb_5kb.mmetacool

metacooler zoomify --binsize 5000 --output assembly_1kb_5kb_2.metacool assembly_1_rep1_DpnII.metacool
metacooler info assembly_1_rep1_DpnII.metacool
metacooler info assembly_1kb_5kb_2.metacool
