library(stringr)

##parsing report function

parsingBowtie2Log<-function(bowtie2log){
  data = readLines(bowtie2log)
  N_reads_mapped=as.numeric(sapply(strsplit(head(data,1), " "), "[[", 1))*as.numeric(str_replace(sapply(strsplit(tail(data,1), " "), "[[", 1),"%",""))/100
  return (N_reads_mapped)
}

parsingBowtie1Log<-function(bowtie1log){
  data = readLines(bowtie1log)
  N_reads_mapped=as.numeric(sapply(strsplit(tail(data,1), " "), "[[", 2))
  return (N_reads_mapped)
}

parsingRemoveDupLog<-function(RdLog){
  data = read.table(RdLog, h=T, sep='\t')
  reads_examined=data$UNPAIRED_READS_EXAMINED + data$READ_PAIRS_EXAMINED
  #reads_examined
  reads_dup= data$UNPAIRED_READ_DUPLICATES + data$READ_PAIR_DUPLICATES  + data$READ_PAIR_OPTICAL_DUPLICATES
  #reads_dup
  N_reads_mapped= reads_examined - reads_dup
  return (N_reads_mapped)
}

