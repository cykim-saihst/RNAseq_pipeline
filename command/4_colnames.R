#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
name = args[1]
file_path = paste0('work/',name,'/',name,'_count.csv')
final_file = paste0('output/',name,'_counts.csv')

count = read.csv(file=file_path,col.names =c('gene',name))

write.csv(count,file=final_file,row.names=FALSE)
