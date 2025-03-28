#!/bin/Rscript

library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]

input <- read.table(input_file, header=FALSE)
names = input[,1];
path = paste0("output/",input[,1],"_counts.csv");
path = as.data.frame(path);
first = read.csv(path[1,],header=TRUE)

for ( i in 2:nrow(path)) {
    file = read.csv(path[i,],header=TRUE)
    temp = inner_join(first, file, by="gene")
    first = temp
    }

write.csv(first,file="output/raw_counts.csv", row.names = FALSE)
