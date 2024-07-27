#!/bin/bash

# Generate Table Files. One for each read.
# Tabex file.ktab LIST is used to view counts

WD=../../output
DIR="$WD"/splitFastqs/HG002_20bp

K="31"

count_kmers() {
    local index="$1"
    mkdir -p ./tmp/"$index"
    for chrom in $(ls $DIR); do
        for region in $(ls "$DIR"/"$chrom"); do
            region_num=$(echo $region | cut -d'_' -f1 | tr -dc '0-9')
            if (( region_num == index )); then
                out_dir="$WD"/kmerTables_20bp/"$K"mers/"$chrom"/"$region"
                mkdir -p $out_dir
                for read in $(ls "$DIR"/"$chrom"/"$region"); do
                    read_out="${read%.*}"
                    FastK -t -k"$K" -T1 -Ptmp/"$index" -N"$out_dir"/"$read_out" "$DIR"/"$chrom"/"$region"/"$read"
                done
            fi
        done
        echo "Completed" $chrom
    done
}

# Bootstrap Multithreading
for i in {0..22}; do 
    count_kmers $i &
done
