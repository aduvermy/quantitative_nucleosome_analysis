#!/usr/bin/env Rscript

library(rtracklayer)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

# test if cmd line if ok : if not, return an error
if (length(args) != 3) {
  stop("3 arguments must be supplied (1 input files, 1 output", call.=FALSE)
} else if (length(args)== 3) {

  bw_IP=args[1]

  ##output
  bw_IP_norm=args[2]

  #norm factor
  norm_factor=as.double(readLines(args[3], n = 1))
}

##import data
bw<-import.bw(bw_IP)


##build occupancy normalized
bw$score<-bw$score*norm_factor


##export bw normalized
export.bw(bw, bw_IP_norm)
