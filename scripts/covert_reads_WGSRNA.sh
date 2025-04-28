#!/bin/bash

# Check if correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file=$1
output_file=$2
bin_size=50000

awk -v bin_size=$bin_size '
{
    chr = $1
    pos = $2
    hp = $3
    read_count = $5
    
    # Calculate bin start and end
    bin_start = int(pos/bin_size)*bin_size + 1
    bin_end = bin_start + bin_size - 1
    key = chr"\t"bin_start"\t"bin_end

    # Sum read counts by haplotype
    if (hp == 1) {
        hp1[key] += read_count
    } else if (hp == 2) {
        hp2[key] += read_count
    }

    # Track minimum SNP position per bin
    if (!(key in min_pos) || pos < min_pos[key]) {
        min_pos[key] = pos
    }
}
END {
    # Print all bins that have either hp1 or hp2 counts
    for (k in min_pos) {
        print k, min_pos[k], (k in hp1 ? hp1[k] : 0), (k in hp2 ? hp2[k] : 0)
    }
}' OFS="\t" "$input_file" | sort -k1,1 -k2,2n > "$output_file"

echo "Processing complete. Output saved to $output_file"
