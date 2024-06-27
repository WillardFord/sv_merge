#!/bin/bash

# Generate Table Files. One for each read.
# tabex file.ktab LIST is used to view counts

WD=/Users/wford/Documents/sv_merge/output
DIR="$WD"/splitFastqs/HG002

K="20"

# I bootstrapped multithreading by just setting i to 0..3 and starting the process in a different terminal window
i=3
mkdir -p ./tmp/"$i"
for chrom in $(ls $DIR); do
    chrom_num=$(echo $chrom | tr -dc '0-9')
    if (( chrom_num % 4 == i )); then
        for region in $(ls "$DIR"/"$chrom"); do
            out_dir="$WD"/kmerTables/"$K"mers/"$chrom"/"$region"
            mkdir -p $out_dir
            for read in $(ls "$DIR"/"$chrom"/"$region"); do
                read_out="${read%.*}"
                fastk -t -k"$K" -T1 -P/Users/wford/Documents/sv_merge/scripts/readFiltering/tmp/"$i" -N"$out_dir"/"$read_out" "$DIR"/"$chrom"/"$region"/"$read"
            done
        done
        echo "Completed" $chrom
    fi
done
