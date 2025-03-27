#!/bin/bash

name=$1
thread=$2

cd work

# Sort BAM by read name
samtools sort -@ ${thread} -n -o $name/${name}_sorted.bam -O BAM $name/${name}_raw.bam

# Filter reads: MAPQ ≥ 20, remove unmapped
samtools view -@ ${thread} -h -q 20 -F 0x4 $name/${name}_sorted.bam -o $name/${name}_sorted_Qcut_mapped.bam

# Fix mate information
samtools fixmate -@ ${thread} -m $name/${name}_sorted_Qcut_mapped.bam $name/${name}_sorted_Qcut_mapped_fixmate.bam

# Sort BAM by coordinate
samtools sort -@ ${thread} $name/${name}_sorted_Qcut_mapped_fixmate.bam -o $name/${name}_sorted_Qcut_mapped_fixmate_sorted.bam

# Remove duplicates
samtools markdup -@ ${thread} -r -s $name/${name}_sorted_Qcut_mapped_fixmate_sorted.bam $name/${name}_final.bam

cd ..

# Run quality control
fastqc -t ${thread} -o work/$name/QC work/${name}/*.bam
multiqc work/$name/QC/* -o work/$name/QC
