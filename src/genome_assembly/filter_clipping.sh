#!/bin/bash

# Arguments:
# $1: Input BAM file
# $2: Output BAM file
# $3: Option to exclude clipped reads ("true" or "false")

# If "true", exclude reads with soft or hard clipping (CIGAR operations H or S)
if [ "$3" == 1 ]; then
    samtools view -h $1 | awk '$6 !~ /H|S/ || $1 ~ /^@/ {print}' | samtools view -Sb - > $2
else
    samtools view -Sb $1 > $2
fi
