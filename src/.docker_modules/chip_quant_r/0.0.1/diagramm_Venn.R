#!/usr/bin/env Rscript

library(ggplot2)
library(ggforce)

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

# test if cmd line if ok : if not, return an error
if (length(args) != 5) {
  stop("5 arguments must be supplied (3 input log files, bowtie version used {1,2} and output)", call.=FALSE)
} else if (length(args)== 5) {
  exclu2fasta=args[1]
  undecided=args[2]
  exclu2fastaCalib=args[3]
  versionBowtie=args[4]
  output=args[5]
}


##parsing function
source("/usr/bin/parsingLog.R")


if (versionBowtie==1)
{
  logPath = undecided
  N_undecided <- parsingBowtie1Log(logPath)
  
  logPath=exclu2fasta
  N_exclusiv2fa <- parsingBowtie1Log(logPath)
  
  logPath=exclu2fastaCalib
  N_exclusiv2faCalib <- parsingBowtie1Log(logPath)
  
}

if (versionBowtie==2)
{
  logPath = undecided
  N_undecided <- parsingBowtie2Log(logPath)
  
  logPath=exclu2fasta
  N_exclusiv2fa <- parsingBowtie2Log(logPath)
  
  logPath=exclu2fastaCalib
  N_exclusiv2faCalib <- parsingBowtie2Log(logPath)
  
}


TT = N_exclusiv2fa + N_exclusiv2faCalib + N_undecided  ## total reads mappable
## % reads from FA
fromFA<-signif((N_exclusiv2fa)/TT*100 , digits = 3)
## % reads from FA_calib
fromFAcalib<-signif((N_exclusiv2faCalib)/TT*100 , digits = 3)
## % reads undecided
fromUndecided<- signif(N_undecided/TT*100,digits = 3)



##data.frame to build venn diagramm
df.venn <- data.frame(x = c(-0.666, 0.666),
                      y = c(-0.5, -0.5),
                      reads = c('from fasta (%)', 'from fasta calibration (%)'))
df.venn.annot <- data.frame(x = c(0, -0.800, 0.800),
                            y = c(-0.5, -0.5, -0.5),
                            prct = c(fromUndecided,fromFA,fromFAcalib))




#Make the plot
svg(filename=output,width=5, height=4)
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1, fill = reads)) +
  geom_circle(alpha = .3, size = 0.5, colour = 'black') +
  coord_fixed() +
  theme_void()+
  theme(legend.position = 'top') +
  annotate("text",x = df.venn.annot$x, y = df.venn.annot$y, label = df.venn.annot$prct, size = 5)
dev.off()

