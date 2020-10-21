#!/usr/bin/env Rscript


##install_packages for docker image


options(repos=structure(c(CRAN="http://cloud.r-project.org")))
install.packages(c("nlme", "Matrix"))
install.packages("https://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-30.tar.gz", dep=T)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rtracklayer")
BiocManager::install("BRGenomics")
install.packages("ggplot2")
install.packages("ggforce")
install.packages("stringr")
