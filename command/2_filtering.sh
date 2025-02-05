#!/bin/bash

name=$1
thread=$2

cd work

# Sort bam by name
samtools sort -@ ${thread} -n -o $name/${name}_sorted.bam -O BAM $name/${name}_raw.bam

# Filter and map quality cut
samtools view -@ ${thread} -h -q 30 -F 0x4 -F 0x400 $name/${name}_sorted.bam -o $name/${name}_sorted_Qcut_mapped.bam

# Fix mate information
samtools fixmate -@ ${thread} -m $name/${name}_sorted_Qcut_mapped.bam $name/${name}_sorted_Qcut_mapped_fixmate.bam

# Sort bam by coordinate
samtools sort -@ ${thread} $name/${name}_sorted_Qcut_mapped_fixmate.bam -o $name/${name}_sorted_Qcut_mapped_fixmate_sorted.bam

# Remove duplicates
samtools markdup -@ ${thread} -r -s $name/${name}_sorted_Qcut_mapped_fixmate_sorted.bam $name/${name}_final.bam

cd ..

fastqc -t ${thread} -o work/$name/QC  work/${name}/*.bam
multiqc work/$name/QC/* -o work/$name/QC
