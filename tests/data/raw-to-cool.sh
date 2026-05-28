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

## Fix all timestamps backed into HDF5 attr `creation-date` to 2001-01-01T00:00:00Z (source date epoch) for reproducibility, and zero out the 4-byte checksum field in pairs.gz files to avoid non-determinism from gzip's default timestamp-based checksum.
find . -name "*.cool" -o -name "*.metacool" -o -name "*.mcool" -o -name "*.mmetacool" -maxdepth 1 | while read f; do
    echo "Fix timestamp in $f"
    python -c "import h5py; import os; from datetime import datetime, timezone; import h5py; f = h5py.File('$f', 'r+'); f.attrs['creation-date'] = datetime(2001, 1, 1, tzinfo=timezone.utc).isoformat(); f.close()"
    micromamba run -n hdf5 h5repack -i "$f" -o "$f.tmp" && mv "$f.tmp" "$f"
done
find . -name "*.mcool" -maxdepth 1 | while read f; do
    echo "Fix timestamp in $f"
    python -c "import h5py; from datetime import datetime, timezone; dt = datetime(2001, 1, 1, tzinfo=timezone.utc).isoformat(); f = h5py.File('$f', 'r+'); [f[f'/resolutions/{res}'].attrs.update({'creation-date': dt}) for res in f['/resolutions'].keys()]; f.close()"    done
    micromamba run -n hdf5 h5repack -i "$f" -o "$f.tmp" && mv "$f.tmp" "$f"
done
find . -name "*.mmetacool" -maxdepth 1 | while read f; do
    echo "Fix timestamp in $f"
    python -c "import h5py; from datetime import datetime, timezone; dt = datetime(2001, 1, 1, tzinfo=timezone.utc).isoformat(); f = h5py.File('$f', 'r+'); [f[f'/entries/{res}'].attrs.update({'creation-date': dt}) for res in f['/entries'].keys()]; f.close()"    done
    micromamba run -n hdf5 h5repack -i "$f" -o "$f.tmp" && mv "$f.tmp" "$f"
done
find . -name "*.pairs.gz" -maxdepth 1 | while read f; do
    echo "Fix checksum in $f"
    printf '\x00\x00\x00\x00' | dd of="$f" bs=1 seek=4 count=4 conv=notrunc
done
