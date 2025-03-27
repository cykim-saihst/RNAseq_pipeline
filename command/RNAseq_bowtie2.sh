#!/bin/bash

name=$1
thread=$2
ref=$3

bash command/0_mkdir.sh $name 
bash command/1_bowtie2.sh $name $thread $ref
bash command/2_filtering.sh $name $thread 
bash command/3_featureCounts.sh $name $thread $ref
Rscript command/4_colnames.R $name 
