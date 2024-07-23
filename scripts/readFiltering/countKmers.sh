#!/bin/bash

# Generate Table Files. One for each read.
# Tabex file.ktab LIST is used to view counts

WD=../../output
DIR="$WD"/splitFastqs/HG002_20bp

K="21"

count_kmers() {
    local index="$1"
    mkdir -p ./tmp/"$index"
    for chrom in $(ls $DIR); do
        chrom_num=$(echo $chrom | tr -dc '0-9')
        if (( chrom_num == index )); then
            for region in $(ls "$DIR"/"$chrom"); do
                out_dir="$WD"/kmerTables_20bp/"$K"mers/"$chrom"/"$region"
                mkdir -p $out_dir
                for read in $(ls "$DIR"/"$chrom"/"$region"); do
                    read_out="${read%.*}"
                    #echo /Users/wford/Documents/sv_merge/scripts/readFiltering/tmp/"$index" "$out_dir"/"$read_out" "$DIR"/"$chrom"/"$region"/"$read"
                    FastK -t -k"$K" -T1 -Ptmp/"$index" -N"$out_dir"/"$read_out" "$DIR"/"$chrom"/"$region"/"$read"
                done
            done
            echo "Completed" $chrom
        fi
    done
}

# Bootstrap Multithreading
for i in {0..22}; do 
    count_kmers $i &
done
