#!/bin/bash

# Check if correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <BAM_FILE> <OUTPUT_FILE>"
    exit 1
fi

# Define input BAM file and output file
BAM_FILE=$1
BINSIZE=$2
OUTPUT_FILE=$3

# Check if BAM file exists
if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file '$BAM_FILE' not found!"
    exit 1
fi

# Process BAM file and count methylation levels per CpG site
samtools view -F 3852 -q 20 "$BAM_FILE" | awk -v output_file="$OUTPUT_FILE" '
BEGIN {
    OFS="\t";
}

{
    chrom = $3;      # Chromosome
    start = $4;      # Start position of read
    seq = $10;       # Read sequence
    has_hp = 0;
    has_ps = 0;
    ps_tag = "NA";
    hp_tag = "NA";
    methylation = "NA";

    # Scan auxiliary fields for HP, PS, and XM tags
    for (i=11; i<=NF; i++) {
        if ($i ~ /^HP:i:/) { has_hp = 1; hp_tag = substr($i, 6); }
        if ($i ~ /^PS:i:/) { has_ps = 1; ps_tag = substr($i, 6); }
        if ($i ~ /^XM:Z:/) { methylation = substr($i, 6); }
    }


    # Process only reads that have HP, PS, and methylation tags
    if (has_hp && has_ps && methylation != "NA") {
        for (j=1; j<=length(methylation); j++) {
            meth_status = substr(methylation, j, 1);
            pos = start + j - 1;  # Genomic position based on start position

            if (meth_status != "Z" && meth_status != "z") continue;

            # Count Z (methylated) and z (unmethylated)
            methylated_count = (meth_status == "Z") ? 1 : 0;
            unmethylated_count = (meth_status == "z") ? 1 : 0;

            # Output: chrom, PS, HP, CpG position, methylated count, unmethylated count
                print chrom, ps_tag, hp_tag, pos, methylated_count, unmethylated_count >> output_file;
        }

    }
}

END {
    close(output_file);
}' 

echo "Methylation counts saved to $OUTPUT_FILE"

