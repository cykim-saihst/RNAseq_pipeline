#!/bin/bash

name=$1
thread=$2
ref=$3

if [ ! -f "ref/${ref}.gtf.gz" ]; then
  echo "❌ Error: Reference annotation file 'ref/${ref}.gtf.gz' not found!"
  exit 1
fi

if [ ! -f "work/$name/${name}_final.bam" ]; then
  echo "❌ Error: BAM file 'work/$name/${name}_final.bam' not found!"
  exit 1
fi


# featureCounts 실행
featureCounts -p --countReadPairs -B -C -a ref/${ref}.gtf.gz -o work/$name/${name}_count work/$name/${name}_final.bam

# 유전자 ID 패턴 결정 (모든 버전의 마우스/사람 지원)
if [[ "$ref" == mm* ]]; then
  gene_prefix="ENSMUSG"
elif [[ "$ref" == hg* ]]; then
  gene_prefix="ENSG"
else
  echo "❌ Error: Unsupported reference genome '$ref'."
  echo "   This pipeline only supports mouse (mmXX) and human (hgXX) genomes."
  exit 1
fi

# 필터링 과정
echo "🔍 Filtering featureCounts output for gene prefix: ${gene_prefix}"
grep "${gene_prefix}" work/$name/${name}_count > temp
awk '{print $1",",$7}' temp > work/$name/${name}_count.csv

rm temp

