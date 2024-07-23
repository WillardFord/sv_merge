#!/bin/bash

# Copy all fastq files to zeppelin

DIR=../../output/HG002_20bp/

split_fastq() {
    local index="$1"
    for chrom in $(ls $DIR); do
        chrom_num=$(echo $chrom | tr -dc '0-9')
        if (( chrom_num % 4 == i )); then
            scp -r "$DIR""$chrom" wford@zep-gpu.broadinstitute.org:/data/wford/sv_merge/output/HG002_20bp/
            echo "Completed" $chrom
        fi
    done
}

# Bootstrap Multithreading
for i in {0..3}; do 
    split_fastq $i &
done
