#!/bin/bash

gene_list="ecm33 SPNCRNA.1452 sol1 rds1 exg1 gas1 grx4 SPNCRNA.1482 adg2 eng1 prl53 rpc82 cdc22 fba1 bos1 SPBC1921.05 tef103 SPNCRNA.229 SPCC1919.05 SPAC26H5.10c pma1 SPBC1E8.05 cnt5 gas5 SPNCRNA.1459 dlc2 htb1 SPAC4H3.12c sme2 cdt1 utp15 SPNCRNA.659 slp1 SPSNRNA.01 mid2 sad1 SPNCRNA.163 SPAC222.08c sfp1 ace2 srp2 taf9 cdc18 cfh4 seb1 SPBC713.10 SPBC713.11c SPNCRNA.906 SPCC1393.08 tos4"

## check if output file already exist. If it does, the file is removed
if [ -e "condensin_peaks_strongest_gene.gff3" ]
then
  rm condensin_peaks_strongest_gene.gff3
fi

## create the output file, empty for now (initialization)
touch condensin_peaks_strongest_gene.gff3


## for each gene, we add the corresponding gff3 line to the output (only the line for gene, not for exon, transcript, ... description)
for xxx in $gene_list
do
  echo $xxx
  res=$(grep $xxx S_pombe_ASM294v226.gff3 | grep -P "PomBase\tgene\t")
  if  [[ -z "${res}" ]]
  then
    res=$(grep $xxx S_pombe_ASM294v226.gff3 | grep -P "_gene\t")
    grep $xxx S_pombe_ASM294v226.gff3 | grep -P "_gene\t" >> condensin_peaks_strongest_gene.gff3
  else
    grep $xxx S_pombe_ASM294v226.gff3 | grep -P "PomBase\tgene\t" >> condensin_peaks_strongest_gene.gff3
  fi

  if [[ -z "${res}" ]]
  then
     echo "nothing"
  else
     echo "$res"
  fi
  
  #echo "$res" >> condensin_peaks_strongest_gene.gff3
done
