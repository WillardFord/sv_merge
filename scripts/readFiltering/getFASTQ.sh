#!/bin/bash

BED_DIR="/Users/wford/Documents/windows_per_sample_per_chrom_flank_0_dipcall_het/intersection_002_733/"

for sample in HG002 HG00733
do
    for i in $(seq 1 22);
    do
        OUTPUTDIR="../output/$sample/chr$i"
        BED="$BED_DIR"chr"$i".bed""
        BAM="gs://fc-28761d6c-5677-4941-86e7-6e42b59a27f4/willard/"$sample"/"$sample"_haplotagged.bam "
        echo $BAM
        echo $BED
        echo $OUTPUTDIR 
        ./extract_reads_from_windows \
            --bam $BAM \
            --bed $BED \
            --output_dir $OUTPUTDIR \
            --flank_length 200 \
            --n_threads 4 \
            --require_spanning \
            --tags PS,HP \
            --force_forward
    done
done