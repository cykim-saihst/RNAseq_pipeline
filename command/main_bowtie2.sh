#!/bin/bash

# Input 파일 요구 사항 안내
echo "========================================="
echo " RNA-seq Preprocessing Pipeline - Bowtie2"
echo "========================================="
echo "This script supports mouse (Mus musculus) and human (Homo sapiens) reference genomes."
echo "Each sample should have two paired FASTQ files in the 'input/' directory."
echo ""
echo "Expected Input File Format: <sample_name> <threads> <reference>"
echo "Example:"
echo "  sample1 5 mm10"
echo "  sample2 3 hg38"
echo ""
echo "All required files must be placed in the following directories:"
echo "  input/ directory:"
echo "    - input/input_file.txt"
echo "    - input/sample1_1.fq.gz"
echo "    - input/sample1_2.fq.gz"
echo "  ref/ directory:"
echo "    - ref/mm10.fa, ref/mm10.gtf.gz"
echo "    - ref/hg38.fa, ref/hg38.gtf.gz"
echo "  command/ directory :"
echo "    - command/0_mkdir.sh"
echo "    - command/1_bowtie2.sh"
echo "    - command/2_filtering.sh"
echo "    - command/3_featureCounts.sh"
echo "    - command/4_colnames.R"
echo "    - command/RNAseq_bowtie2.sh"
echo "    - command/main_bowtie2.sh"
echo "    - command/merge_output.R"
echo "-----------------------------------------"

# 사용자에게 input file 입력 요구
while true; do
  read -p "Enter the input file name (located in 'input/' directory): " input_file
  input_file="input/${input_file}"  

  if [ -f "$input_file" ]; then
    break
  else
    echo "❌ Error: File '$input_file' not found in 'input/' directory! Please enter a valid file name."
  fi
done

# Input 파일을 읽으면서 존재 여부 체크 및 실행
echo "Checking required input files..."
while read -r input1 input2 input3 || [ -n "$input1" ] || [ -n "$input2" ] || [ -n "$input3" ]; do
  file1="input/${input1}_1.fq.gz"
  file2="input/${input1}_2.fq.gz"
  ref_fa="ref/${input3}.fa"
  ref_gtf="ref/${input3}.gtf.gz"

  # 유전체가 마우스인지 사람인지 판별
  if [[ "$input3" == mm* ]]; then
    species="mouse"
  elif [[ "$input3" == hg* ]]; then
    species="human"
  else
    echo "❌ Error: Unsupported reference genome '$input3'."
    echo "   This pipeline only supports mouse (mmXX) and human (hgXX) genomes."
    echo "-----------------------------------------"
    exit 1
  fi

  # 파일 존재 여부 확인
  if [ ! -f "$file1" ] || [ ! -f "$file2" ] || [ ! -f "$ref_fa" ] || [ ! -f "$ref_gtf" ]; then
    echo "❌ Error: Required files not found for sample '$input1'."
    echo "  Missing files:"
    [ ! -f "$file1" ] && echo "  - $file1"
    [ ! -f "$file2" ] && echo "  - $file2"
    [ ! -f "$ref_fa" ] && echo "  - $ref_fa (Reference genome FASTA file not found)"
    [ ! -f "$ref_gtf" ] && echo "  - $ref_gtf (Reference annotation file not found)"
    echo "-----------------------------------------"
    exit 1
  fi
done < "$input_file"

echo "✅ All required input files are present. Starting analysis..."
echo "========================================="

# 본격적으로 분석 실행
while read -r input1 input2 input3 || [ -n "$input1" ] || [ -n "$input2" ] || [ -n "$input3" ]; do
  echo "Processing: input/${input1}_1.fq.gz and input/${input1}_2.fq.gz with ${input2} threads (Reference: ${input3})"
  bash command/RNAseq_bowtie2.sh "$input1" "$input2" "$input3"
done < "$input_file"

Rscript command/merge_output.R "$input_file"
