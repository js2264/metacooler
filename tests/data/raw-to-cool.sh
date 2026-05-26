#! /usr/bin/env bash

## Create cools
for i in 1 2; do
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 assembly_${i}.chrom.sizes:1000 contacts_${i}_rep1.pairs.gz assembly_${i}_rep1_1kb.cool
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 assembly_${i}_DpnII-digested.bed contacts_${i}_rep1.pairs.gz assembly_${i}_rep1_DpnII.cool
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 assembly_${i}.chrom.sizes:1000 contacts_${i}_rep2.pairs.gz assembly_${i}_rep2_1kb.cool
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 assembly_${i}_DpnII-digested.bed contacts_${i}_rep2.pairs.gz assembly_${i}_rep2_DpnII.cool
    cooler zoomify --balance -r 1000,2000,5000 assembly_${i}_rep1_1kb.cool
    cooler zoomify --balance -r 1000,2000,5000 assembly_${i}_rep1_DpnII.cool
done

## Create metacool from pairs
for i in 1 2; do
    metacooler create pairs --force --contigs contigs_table_${i}.tsv --genomes mag_table_${i}.tsv --pairs contacts_${i}_rep1.pairs.gz --bins 1000 --output assembly_${i}_rep1_1kb.metacool
    metacooler create pairs --force --contigs contigs_table_${i}.tsv --genomes mag_table_${i}.tsv --pairs contacts_${i}_rep2.pairs.gz --bins 1000 --output assembly_${i}_rep2_1kb.metacool
    metacooler create pairs --force --contigs contigs_table_${i}.tsv --genomes mag_table_${i}.tsv --pairs contacts_${i}_rep1.pairs.gz --bins assembly_${i}_DpnII-digested.bed --output assembly_${i}_rep1_DpnII.metacool
    metacooler create pairs --force --contigs contigs_table_${i}.tsv --genomes mag_table_${i}.tsv --pairs contacts_${i}_rep2.pairs.gz --bins assembly_${i}_DpnII-digested.bed --output assembly_${i}_rep2_DpnII.metacool
done

## Create mmetacool
metacooler merge --force --output assembly_1_1kb.metacool assembly_1_rep1_1kb.metacool assembly_1_rep2_1kb.metacool
metacooler stack --force --output assembly_1_1kb.mmetacool assembly_1_rep1_1kb.metacool assembly_1_rep2_1kb.metacool
