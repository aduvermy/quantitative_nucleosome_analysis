#!/usr/bin/env Rscript


library(rtracklayer)
library(BRGenomics)
library(ggplot2)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

# test if cmd line if ok : if not, return an error
if (length(args) != 3) {
  stop("3 arguments must be supplied (2 input files and 1 output).n", call.=FALSE)
} else if (length(args)== 3) {
  ##log mapping
  bw_IP_C_wt=args[1]
  bw_IP_C_mut=args[2]
  ##output
  output=args[3]
}

## import BigWIg and merged
IP_C_mut<-import.bw(bw_IP_C_mut)
IP_C_mut<-makeGRangesBRG(IP_C_mut,ncores = 1)

IP_C_wt<-import.bw(bw_IP_C_wt)
IP_C_wt<-makeGRangesBRG(IP_C_wt,ncores = 1)

dta_C <- mergeGRangesData(IP_C_wt, IP_C_mut, multiplex = TRUE,ncores = 1)
rm(IP_C_wt, IP_C_mut) ##rm useless variable

## sub sample
dta_C <- as.data.frame(dta_C)
dta_C <- dta_C[sample(nrow(dta_C), 20000), ]

## calcul mean ratio
meanRatio=mean(dta_C$IP_C_wt)/mean(dta_C$IP_C_mut)

##build plot
svg(filename=output, width=10, height=10)
ggplot(data= dta_C, aes(x = IP_C_wt/IP_C_mut, y=stat(density)))+
  geom_histogram(color="black",fill = "white", bins=50)+
  geom_density(alpha = 0,size=1)+
  geom_vline(aes(xintercept = meanRatio, col=toString(meanRatio)), linetype = "dashed", size = 1)+
  labs(col = "mean(IP_C_wt)/mean(IP_C_mut)")
dev.off()
