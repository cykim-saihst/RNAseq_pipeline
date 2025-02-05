#!/bin/bash

name=$1
thread=$2
ref=$3

# Bowtie2 index 확인 후 없으면 생성
index_prefix="ref/${ref}"
if [ ! -f "${index_prefix}.1.bt2" ]; then
  echo "🔍 Bowtie2 index not found for reference '${ref}', creating index..."
  bowtie2-build ref/${ref}.fa ref/${ref}
  echo "✅ Bowtie2 index created and stored in 'ref/${ref}.*.bt2'"
else
  echo "✅ Bowtie2 index already exists for '${ref}'. Skipping index generation."
fi

# Bowtie2 실행
echo "🔬 Running Bowtie2 for sample '${name}'..."
bowtie2 -p ${thread} -x ref/${ref} --sensitive -X 2000 -1 input/${name}_1.fq.gz -2 input/${name}_2.fq.gz > work/$name/${name}_raw.sam 

# SAM → BAM 변환
samtools view -@ ${thread} -b -S -o work/$name/${name}_raw.bam work/$name/${name}_raw.sam 

# SAM 파일 삭제
rm work/$name/${name}_raw.sam
