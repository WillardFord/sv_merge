#!/bin/bash

# Split a fastq file into reads with haplotypes
split_fastq() {
    local input_fastq="$1"
    local output_directory=$(get_output $input_fastq)

    mkdir -p "$output_directory"
    declare -a reads_with_numbers
    # Read 4 lines at a time from the input FASTQ file
    while IFS= read -r line1 && IFS= read -r line2 && IFS= read -r line3 && IFS= read -r line4; do
        # Only take lines that have assigned haplotypes
        if [[ $line1 =~ HP:i:([0-9]+) ]]; then
            hp_number="${BASH_REMATCH[1]}"  # Extract the number after HP:i:
            reads_with_numbers+=("$hp_number|$line1|$line2|$line3|$line4")
        fi
    done < "$input_fastq"

    # Perform stable sort based on the extracted number
    IFS=$'\n' sorted_reads=($(sort -t '|' -k1 -s -n <<<"${reads_with_numbers[*]}"))

    local count=0
    for read_info in "${sorted_reads[@]}"; do
        ((count++))
        # Extract lines from sorted read_info
        IFS='|' read -r hp_number line1 line2 line3 line4 <<< "$read_info"
        {
            echo "$line1"
            echo "$line2"
            echo "$line3"
            echo "$line4"
        } > "$output_directory/read$count.fastq"
    done
}

# Get output directory
get_output() {
    local input_fastq="$1"
    local file_name="$(basename $input_fastq)"
    local chr_name="$(basename $(dirname $input_fastq))"
    local sample_name="$(basename $(dirname $(dirname $input_fastq)))"
    
    local suffix="$sample_name"/"$chr_name"/"$file_name"
    local prefix=${input_fastq:0:12}

    local output_directory="$prefix"/splitFastqs/"$suffix"
    output_directory="${output_directory%.*}"
    echo $output_directory
}

sample_dir=../../output/HG002_20bp
for chrom_dir in $(ls $sample_dir); do
    for fastq in $(ls "$sample_dir"/"$chrom_dir"); do
        input_fastq="$sample_dir"/"$chrom_dir"/"$fastq"

        split_fastq "$input_fastq" &
    done
    echo "Completed" $chrom_dir
done
