#!/usr/bin/env Rscript

#library(data.table)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

# test if cmd line if ok : if not, return an error
if (length(args) != 5) {
  stop("5 arguments must be supplied (4 input files, and tools used {1,2}).n", call.=FALSE)
} else if (length(args)== 5) {
  ##log mapping
  IP_X=args[1]
  IP_C=args[2]
  WCE_X=args[3]
  WCE_C=args[4]
  ##tools used for log
  tools=args[5]

}



##functions for parsing Log reports
source("/usr/bin/parsingLog.R")

if (tools==1)
{
  logPath = IP_X
  N_IP_X <- as.double(parsingBowtie1Log(logPath))

  logPath = IP_C
  N_IP_C <- as.double(parsingBowtie1Log(logPath))

  logPath = WCE_X
  N_WCE_X <- as.double(parsingBowtie1Log(logPath))

  logPath = WCE_C
  N_WCE_C <- as.double(parsingBowtie1Log(logPath))
}

if (tools==2)
{
  logPath = IP_X
  N_IP_X <- as.double(parsingBowtie2Log(logPath))

  logPath = IP_C
  N_IP_C <- as.double(parsingBowtie2Log(logPath))

  logPath = WCE_X
  N_WCE_X <- as.double(parsingBowtie2Log(logPath))

  logPath = WCE_C
  N_WCE_C <- as.double(parsingBowtie2Log(logPath))
}


if (tools=="RD")
{
  logPath = IP_X
  N_IP_X <- as.double(parsingRemoveDupLog(logPath))

  logPath = IP_C
  N_IP_C <- as.double(parsingRemoveDupLog(logPath))

  logPath = WCE_X
  N_WCE_X <- as.double(parsingRemoveDupLog(logPath))
  
  logPath = WCE_C
  N_WCE_C <- as.double(parsingRemoveDupLog(logPath))
}


OR<-N_WCE_C*N_IP_X/(N_WCE_X*N_IP_C)
cat(OR)


