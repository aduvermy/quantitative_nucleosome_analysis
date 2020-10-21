# STAGE ENS

Compte rendu avancé projet


## Resumé biblio/Intro

  Chez la plupart des eucaryotes, les fibres de chromatine se métamorphosent et se compactent pour former des chromosomes individuels, pendant la mitose et la meiose.
Cette profonde réorganisation, appelée condensation des chromosome, est nécessaire pour la bonne ségrégation des chromosomes. De la levure aux mammifères, cette condensation
implique les topoisomérases mais aussi le complexe condensine. Selon certaines études, il semble que condensine modifie la topologie de l'ADN en introduisant des superenroulement en attrapant
la molecule d'ADN dans sa structure en anneau.

  Au cours du cycle cellulaire, les nucleosomes jouent aussi un rôle important dans l'organisation du génome eucaryote, en facilitant son repliement au sein du noyau.
L'organisation des nucléosomes module l'accès à l'ADN et régulent ainsi des processus biologiques basiques tels que la transcription, la réplication
ou la recombinaison. Par séquençage MNase, les travaux de Soriano et al, ont pu montrer, sous différentes conditions physiologiques, que le profil des nucleosomes de S.pombe est largement
maintenu sur tout le génome, et notamment dans les régions déplétées en nucléosomes (NDRs) présentent immédiatement en amont du début de la transcription de la plupart des gènes.
Ils ont ainsi pu cartographier, le long du génome, les NDRs de S.pombe et montrer que les NDRs sont des clusters de liaison entre l'ADN et les FTs.

  Les travaux de Toselli et al, ont rapporté que condensine s'accumulait à proximité des régions déplétés en nucléosomes(NDRs), durant la mitose des levures. Ils ont aussi proposé le fait que
les nucléosomes agissent comme une barrière pour la liaison de condensine et par conséquent que les NDRs formés au niveau des régions hautement traduites constitue des points d'accès d'ADN libre pour condensine.

  Afin de corroborer cette hypothèse nous utiliserons l'outils DANPOS, les données MNase de Soriano et al, pour cartographier les NDRs le long du génome de S.pombe sur le modèle
de l'analyse de Soriano et al,.  
Des données ChipSeq condensine de Kakui et al, et de Sutani permettront ensuite d'identifier les zones d'interactions entre l'ADN et condensine. Ces zones d'interactions pourront alors être comparées aux régions déplétées en nucleosomes.  
Plusieurs études ont montré qu’en l’absence de contrôle interne invariable (spike-in), la ChIP-seq ne fournit pas une mesure quantitative de l’association d’un facteur au génome (Bonhoure et al., 2014; Hu et al., 2015; Orlando et al., 2014). Sur la base d’une méthode récemment publiée (Hu et al., 2015), nous avons développé un protocole de ChIP-Seq calibré en utilisant de la chromatine de S. cerevisiae comme référentiel de calibration (spike-in control). En s’inspirant des approches décrites dans les travaux récents (Bonhoure et al., 2014; Hu et al., 2015; Jeronimo et al., 2019; Larochelle et al., 2018; Orlando et al., 2014) nous tenterons de quantifier l'interaction de condensine à l'ADN en utilisant comme spike-in control une ChiP seq de Rpb1 chez cerevisiae publiée par Larochelle et al.

## Mat&Meth

**1) Analyse de données MNase-seq utilisant DANPOS (algo Dpos) version 2.2.2 pour mettre en évidence les NDRs à partir des données de Toselli.**  
**2) Analyse de données ChIP-seq utilisant MACS2 pour mettre en évidence les peaks condensine**  
**2) Comparer les positions des NDRs avec celles de pics condensin en mitose  avec Hawkes**  
**3) Mettre en place une analyse de ChIP-seq quantitatif**

**Pour garantir la reproductibilité et la pérennité de l'analyse le pipeline sera intégré à Nextflow, en utilisant des images docker.**  

Données utilisées:    
- *Toselli-Mollereau et al., 2016* (Mnaseq S.pombe)    

|       SRA        |        TYPE        |        SE OU PE        |
|:----------------:|:------------------:|:----------------------:|
| ERR1353027       | MNASE              |          SE            |
| ERR1353028       | MNASE              |            SE          |
| ERR1353029       | MNASE              |           SE           |


- *Kakui et al., 2017* (ChIP seq condensine S.pombe)  

|       SRA        |        TYPE        |        SE OU PE        |
|:----------------:|:------------------:|:----------------------:|
| SRR5227705       | chip (input)       |          PE            |
| SRR5227706       | chip condensine    |            PE          |

- *Sutani et al.* (ChIP seq condensine S.pombe)  

|       Data set   |   SRA              |         Type           |  INPUT ou IP       | Séquençage     |
|:----------------:|:------------------:|:----------------------:|:------------------:|:--------------:|
| Data set 1  | SRX691442: S_pombe_NLS-GFP-PK_ChIP_input_mitosis | Control de spécificité (CHIP contre une protéine GFP-Pk nucléaire nucléosoluble. Control de spécificité pour les CHIP condensin-Pk | Input | Illumina|
| Data set 1  | SRX691441: S_pombe_NLS-GFP-PK_ChIP_mitosis |   Control de spécificité (CHIP contre une protéine nucléaire nucléosoluble)  | IP | Illumina |
| Data set 2  | SRX687805: no_tag_ strain_ChIP_input_mitosis |  Control de bruit de fond pour les ChIP anti-PK. Une ChIP avec un anti-PK est réalisée sur de la chromatine dépourvue de tag Pk9 | Input | AB solid|
| Data set 2| SRX685977: no_tag_strain_ChIP_mitosis | Control de bruit de fond pour les ChIP anti-PK | IP | AB SOLID |
| Data set 3 |SRX691444: S_pombe_RNApol2_ChIP_input_mitosis |ChIP RNA Pol II| INPUT | AB SOLID |
| Data set 3 |SRX691443: S_pombe_ RNApol2_ChIP_mitosis | ChIP RNA Pol II | IP | AB SOLID |
| Data set 4 |SRX674032: S_pombe_Cut14-PK_ChIP_mitosis_1 | ChIP condensin (Cut14-PK) IP anti-PK | IPreplicat 1 | AB SOLID |
| Data set 4 |SRX674033: S_pombe_Cut14-PK_ChIP_mitosis_1 | ChIP condensin (Cut14-PK) IP anti-PK | INPUT replicat 1 | AB SOLID |
| Data set 5 |SRX674035: S_pombe_Cut14-PK_ChIP_mitosis_2 | ChIP condensin (Cut14-PK) IP anti-PK | IP replicat 2 | AB SOLID |
| Data set 5 |SRX674037: S_pombe_Cut14-PK_ChIP_mitosis_2 | ChIP condensin (Cut14-PK) IP anti-PK | INPUT replicat 2 | AB SOLID |

- *Larochelle et al., 2018* (chip quantit)  

|       SRA        |        TYPE        |        SE OU PE        | S.pombe/S.cerevisiae |
|:----------------:|:------------------:|:----------------------:|:--------------------:|
| SRR7289786       | IP WT              |          SE            | 9:1                  |
| SRR7289787       | WCE WT             |            SE          | 9:1                  |
| SRR7289788       | IP MUT             |           SE           |9:1                   |
| SRR7289789       | WCE MUT            |           SE           |9:1                   |


- Genome Reference S.Pombe  
      - **Ensembl ASM294v2, May 2009 GCA_000002945.2**  

- Genome Reference S.Cerevisiae  
      - **GCA_000146045.2 from ENA/EMBL full**


## Pipeline MNase-seq vs Toselli

|       Step       |   mnase_analysis   |         Toselli        |
|:----------------:|:------------------:|:----------------------:|
| adapter removal  | None               |          None          |
| trimming         | None               |               None     |
| mapping          | Bowtie2            |  Bowtie2               |
| duplicate removal| Oui                |None                    |
| Peak calling     | DANPOS             |Matlab “Wavelet Toolbox”(Soriano)|
| Normalisation    |average occupancy   |average occupancy       |
| NDRs             |L>150, occ_norm>0.4 |L>150, occ_norm>0.4     |


*N.B: Pour mes résultats, j'ai enlevé les reads dupliqués (avec picardTools)*
*Il semble que cela n'ai pas d'effet sur les résultats (même nb de NDRs avec/sans)*

**In order to use our pipeline please visit README_mnase.md**

## Pipeline chIP-seq vs Kakui

|       Step       |  chip_analysis  |       Kakui       |
|:----------------:|:---------------:|:-----------------:|
| adapter removal  | Oui             |Non                |
| trimming         | None            |oui to 50 bp       |
| mapping          | Bowtie2         |bwa mem            |
| duplicate removal| Oui             |None               |
| Peak calling     | MACS2           |MACS2              |

*N.B: Pour mes résultats, j'ai enlevé les reads dupliqués (avec picardTools)*
*Pas testé, mais théoriquement pas d'effet sur les résultats puisque MACS2 ne retient qu'une lecture par position afin d'éviter les bais d'amplifications*

**In order to use our pipeline please visit README_chip.md**


## Pipeline chIP-seq vs Sutani


|       Step       |  chip_analysis  |       Sutani       |
|:----------------:|:---------------:|:------------------:|
| adapter removal  | non             |Non                 |
| trimming         | Non             | Non                |
| mapping          | Bowtie1         |bowtie1             |
| duplicate removal| Oui             |None                |
| Peak calling     | MACS2           |parse2wig and DROMPA|

**N.B:‘Fewer paired peaks X than 1000’ means that MACS only identified X model peaks and may indicate potential data quality issues because 1,000 model peaks are needed to robustly estimate ChIP-DNA fragment size**  
-> utlisation fix-bimodal option to skip building model if needed  
*Macs2 n'arrivait pas à construire le modèle (estimation de la taille des fragments) l'ajout de cette option a donc été nécessaire pour les données Sutani (SE)*

**In order to use our pipeline please visit README_chip.md**


## Analyse quantitative

0)QC  
  -> reads bonne qualité  
1)MAPPING pour déterminer reads spécifique d'une espèce  
2)Calcul du ratio occupation (défini par Hu)  s
3)Normalisation de l'occupation  
  -> Objectif tenté de faire fiter au max les IP Calib mut et wt pour pouvoir comparer position par position l'occupation WT vs MUT des IP X  
4)Comparaison des positions 2 à 2  

## Résultats analyse de données MNAseq


1) mapping

|    ech          |   reads totaux  |    reads mappés   | reads non mappés |
|:---------------:|:---------------:|:-----------------:|:----------------:|
| ERR1353027      | 47 019 988      |46 608 441 (99.12%)|  411547 (0.88%)  |    
| ERR1353028      | 53 548 393      |53 027 827 (99.03%)|  520566 (0.97%)  |
| ERR1353029      | 53 201 667      |52 745 553 (99.14%)|  456114 (0.86%)  |


2) positions des NDRs  

-> Average occupancy normalization

|                 |    Sebastien    |    ERR1353027     |   ERR1353028   |  ERR1353029   | union | intersect |
|:---------------:|:---------------:|:-----------------:|:--------------:|:-------------:|:-----:|:---------:|
| nb NDRs I,II,III| 5 880           |6 667              |  7 672         |   6 031       | 8 368 |  5381     |

*N.B: cmd to union: cat ERR1353027.bed ERR1353028.bed ERR1353029.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > NDRs.bed*  
*N.B: cmd to intersect : bedtools multiinter -i ERR1353027.bed ERR1353028.bed ERR1353029.bed | awk '$4 == 3'  | cut -f1,2,3*  
*/!\ bedtools intersect -a <ref> -b <others files> donne res différents selon l'ordre des fichiers*   
*multiinter permet intersect plusieurs bed et d'éviter effet de reférence. On conserve les régions présentes dans touts les réplicats -> très stringent !*

-> Intersection    

|       Intersect       |   NDRs_union    | NDRs_intersect  |
|:---------------------:|:---------------:|:---------------:|
|     Sebastien         |     5 617       |     4612        |


-> analyse spatiale
![ERR1353027_ERR1353028_ERR1353029_sebastien_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/835e1e850d5ddd6f3f48a0840634cc28/ERR1353027_ERR1353028_ERR1353029_sebastien_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

**L'analyse Hawkes et l'intersection révèlent peu de différences entre les positions des NDRs des 3 réplicats WT.  
Soriano 2013: le profil des nucleosomes de S.pombe est largement maintenu sur tout le génome, et notamment dans les régions déplétées en nucléosomes (NDRs).**  

**Pour la suite 2 approches concernant la liste des NDRs, une plutôt conservative basé sur l'union des 3 réplicats, l'autre plus stringente avec l'intersection.**

## Résultats analyse de données ChIP

0) Reads QUALITY

**N.B: Qualité des données SOLID de Sutani pas fameuse, se ressent sur le taux de mapping (environ 50%)**

![fastqc_per_base_sequence_quality_plot.svg](/uploads/584bc795355837a166b418ce3943ff42/fastqc_per_base_sequence_quality_plot.svg)

1) mapping

|   data               |       ech       | reads ttx |   reads mappés    | reads non mappés  |   Type    |
|:--------------------:|:---------------:|:---------:|:-----------------:|:-----------------:|:---------:|
| Kakui-condensine     |SRR5227705(input)| 62 853 133|53 597 076 (81.76%)|  9256057 (18.24%) |Illumina PE|
| Kakui-condensine     |SRR5227706(chip) |42 676 611 |38 256 277 (89.64%)|  4420334 (10.36%) |Illumina PE|
|Sutani-protNucleosolub|SRR1564297(chip) | 6 880 537 |5 611 914 (81.56%) | 1 268 623 (18.44%)|Illumina SE|
|Sutani-protNucleosolub|SRR1564298(input)| 7 279 549 |7 242 061 (99.49%) |37 488 (0.51%)     |Illumina SE|
|  Sutani-condensine1  |SRR1557175(chip) | 15 093 827|8 437 707 (55.90%) |6 656 120 (44.10%) |SOLID SE   |    
|  Sutani-condensine1  |SRR1557176(input)| 17 870 654|9 427 909 (52.76%) |8 442 745 (47.24%) |SOLID SE   |
|  Sutani-condensine2  |SRR1557178(chip) | 21 626 672|7 189 230 (33.24%) |14 437 442 (66.76%)| SOLID SE  |   
|  Sutani-condensine2  |SRR1559300(input)| 23 191 107|9 590 080 (41.35%) |13 601 027 (58.65%)| SOLID SE  |
|  Sutani-noTag        |SRR1559301(chip) | 8 697 948 |3 601 807 (41.41%) | 5 096 141 (58.59%)|  SOLID SE |
|  Sutani-noTag        |SRR1564296(input)| 12 429 625|9 055 875 (72.86%) |3 373 750 (27.14%) |  SOLID SE |
|  Sutani-RNApolII     |SRR1564402(chip) |22 623 891 |13 112 148 (57.96%)| 9 511 743 (42.04%)| SOLID SE  |
|  Sutani-RNApolII     |SRR1564303(input)|21 961 496 |11 701 166 (53.28%)|10 260 330 (46.72%)|  SOLID SE |
|Larochell-RNApolII-WT |SRR7289786(chip) |21 000 000 |13 000 000 (60%)   | 8 000 000 (40%)   |ILLUMINA SE|    /!\ contient des reads cerevisiae
|Larochell-RNApolII-WT |SRR7289787(input)|8 000 000  |7 800 000 (99%)    | 200 000(1%)       |ILLUMINA SE|    


*Stat mapping Sutani papier: https://static-content.springer.com/esm/art%3A10.1038%2Fncomms8815/MediaObjects/41467_2015_BFncomms8815_MOESM520_ESM.pdf*  
**Nos stats sont légèrement inférieures aux leurs, sans doute lié au version de mapper différent + option plus stringentes.**

2) Peak calling  

*N.B: Toselli-Mollereau: We identified ~7,000 NDRs in mitotic chromosomes in cells arrested in pro/metaphase, but solely ~400 condensin peaks (48 high and 340 low occupancy) have been identified by ChIP-seq at a similar cell cycle stage (Sutani et al, 2015)*
*Fichier comparatif NDRs vs condensineBS.xlsx -> 96 peaks (sebastien_condensine)*

|   data               |       ech       | nb Peak   |
|:--------------------:|:---------------:|:---------:|
| Kakui-condensine     |SRR5227706(chip) | 493       |
|Sutani-protNucleosolub|SRR1564297(chip) | 513       |
|  Sutani-condensine1  |SRR1557175(chip) | 1465      |    
|  Sutani-condensine2  |SRR1557178(chip) | 1123      |   
|  Sutani-noTag        |SRR1559301(chip) | 141       |
|  Sutani-RNApolII     |SRR1564402(chip) |   3167    |
|Larochelle-RNApolII-WT|SRR7289786(chip) |   1897    |
| condensine_sebastien |      (chip)     |   96      |

3) Sutani condensine vs TAG/protNucleosolub

**Les vrais signaux ChIP Cut14-PK (condensine) sont identifiés en comparant [IP Cut14-Pk/Input Cut14-Pk] avec [IP no-tag /Input no-tag] et aussi avec [IP NLS-GFP-Pk/Input NLS-GFP-Pk]. Seuls les pics présents dans l’IP Cut14-Pk et absents dans les deux autres ChIPs contrôles sont validés.**  

*N.B: cmd : bedtools intersect -v -a condensine_rep{1,2}.bed -b noTag.bed protNucleosolub.bed*  

|   data               |       ech       | nb Peak reel   |
|:--------------------:|:---------------:|:--------------:|
|  Sutani-condensine1  |SRR1557175(chip) | 1181           |    
|  Sutani-condensine2  |SRR1557178(chip) | 886            |  

4) Intersect Sutani replicat 1 avec replicat 2

|       Intersect        |Sutani-condensine2_real |
|:----------------------:|:----------------------:|
| Sutani-condensine_real|        693             |


5) Comparaison par intersection

|       Intersect        |Sutani-condensine_real  | condensine_sebastien   |
|:----------------------:|:----------------------:|:----------------------:|
|     kakui              |        247             |       88               |
|   condensine_sebastien |      44                |       96               |



**Comme pour les NDRs, 2 approches concernant les pics condensines une plutôt conservative basée sur l'union de Sutani et Kakui, l'autre plus stringente avec l'intersection.**


|                        |Sutani_intersect_kakui  | Sutani_union_kakui |
|:----------------------:|:----------------------:|:------------------:|
| nb peaks               |        247             |       936          |

*N.B: cmd to intersect: bedtools intersect -a kakui -b Sutani
cmd to union: cat kakui sutani | sort -k1,1 -k2,2n | bedtools merge -i -*


Les données de Sutani se révèle très bruitées.  
Pour les analyses spatiale qui vont suivre nous utilisons uniquement les 493 pics de Kakui.

6) Information complémentaire Kakui/NDRS/polII vs proteine Nucleosoluble

| intersect                 |    kakui(493)    |   RNApolIISutani | RNApolIILarochelle WT| NDRS_UNION(8368)  |NDRS_INTERSECT(5381)|
|:-------------------------:|:----------------:|:----------------:|:--------------------:|:-----------------:|:------------------:|
|proteine Nucleosoluble(513)|        195       |        142       |       123            |       410         |      339           |

-> Analyse spatiale

rouge:condensine vs protNucleosolub  
vert: NDRs vs protNucleosolub  
bleu: RNApol vs protNucleosolub

- Avec NDRS_UNION & RNApolIILarochelle

![NDRs_union_SRR1564297_protNucleosolub_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/e8dc11b05d6d3ec7d00000900d69a74c/NDRs_union_SRR1564297_protNucleosolub_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)


- Avec NDRS_UNION & RNApolIISutani

![NDRs_union_SRR1564297_protNucleosolub_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/d36cd9976d84a735f62ddc53736aad85/NDRs_union_SRR1564297_protNucleosolub_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

195 régions intersectent entre les resultats du chip de la proteine Nucleosoluble et le chip condensine de Kakui, le pics en 0 de l'analyse spatiale conforte cette colocalisation. (pic rouge)
Une 100aine de pics RNApolII (sutani & larochelle) intersectent avec les pics de la prot nucléosoluble, l'analyse spatiale ne fait pas apparaitre de pic en 0.(vallée bleu)
Les pics prot Nucleosoluble intersectent avec les pics RNApolII, mais certainement moins étroitement qu'avec ceux de condensine.  

On peut suspecter que les pics qui intersectent sont le résultat de la fixation aléatoire des proteines étudiées sur l'ADN.(pics incertains).  
**En toute rigueur ces pics devraient être retirés ! (Pas fait ici)**  

Pour ce qui est de l'intersection avec les NDRs, une large majorité des pics "proteine Nucleosoluble" tombent dans les NDRs. Pas étonnant puisque les interactions proteine/ADN ont lieu dans les régions accessibles de l'ADN, soit dans les régions déplétées en nucléosomes.

## Condensine vs NDRs

-> Intersection   

|       Intersect       |       condensine Kakui |      condensine Sutani |  condensine_intersect | condensine_union |
|:---------------------:|:----------------------:|:----------------------:|:---------------------:|:----------------:|
|     NDRs union        |        480/493         |        162/693         |        101            |      541         |
|   NDRs intersection   |        381/493         |         93/693         |      62               |         412      |

## PolII vs NDRs vs condensine

-> Intersection

|       Intersect       |RNApolII Larochell (1897)|
|:---------------------:|:-----------------------:|
|RNApolII Sutani(3167)  |        1509             |

**Les chip-seq RNAPOL II de Sutani et Larochelle concorde du point de vue de leur résultats**  
**-> Majorité des pics de Larochelle retrouvés chez Sutani**  
**Comme pour condensine Sutani RNApol II semble très bruité -> SOLID**


|       Intersect       | condensine Kakui | condensine Sutani |  condensine_intersect |  condensine_union | NDRs intersect | NDRs union |  NDRs Sebastien  |
|:---------------------:|:----------------:|:-----------------:|:---------------------:|:-----------------:|:--------------:|:----------:|:----------------:|
|Sutani-RNApolII (3167) |   335/493        |   619/694         |         221/247       |        733/936    |     1471/5381  |  2385/8368 |  1692/5881       |
|Laroch-RNApolII(1897)  |   227            |   490             |         184           |        583        |     578        |  999       |  636             |

**RNA pol semble largement colocaliser avec condensine.**  
**RNApol semble aussi avoir un lien avec les NDRs, mais pas systématique**  

-> Analyse spatiale avec RNApol Larochelle

rouge: Condensine vs NDRs  
vert: RNApol vs NDRs  
bleu: RNApol vs condensine  

- Avec NDRs union

![NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/81c0b0b060981d7a4af69df7a694886f/NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

- Avec NDRs intersect

![NDRs_intersect_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/fcf0c82d01838eca501ee4ffaea5634e/NDRs_intersect_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

-> Analyse spatiale avec RNApolIISutani

- Avec NDRs union

![NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/b42068591a45e9d15a682840c9834968/NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

- Avec NDRs intersect

![NDRs_intersect_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/aeea628fff76b4f7aa79e6612ea2bbfa/NDRs_intersect_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)


**On recense ~500 peaks condensine pour environ 5 à 8000 NDRS sur tout le génome S.pombe.
Au vue de leurs nombres respectifs de l'intersection et des résultats de la correlation spatiale tous les NDRs ne sont pas associés à un peak de condensine.
En revanche l'analyse spatiale de Hawkes et l'intersection révèlent que les positions des  peaks de condensine sont très fortement corrélés aux positions des NDRs. Ce qui nous conforte dans l'idée que condensine se lie à l'ADN dans les régions déplétés en nucléosomes.**

**L'analyse chip révèle ~3000 pics pour la RNApolII. L'intersection des pics condensine (ttes approches) avec RNApolII, montre une forte colocalisation entre les 2 protéines.  
Lors de ses analyses Sutani a conclu que  RNApolII et condensine colocalisaient. Nos résultats d'intersection corroborent cette conclusion. Mais pas l'analyse hawkes qui révèle une vallée en 0 lors de la comparaison condensine vs polII, finalement contraire aux observations de l'Intersection.  
On notera que les pics condensine et RNApolII sont très larges. Ainsi on peut dire que  les régions associées aux pics obtenus en sortie du pipeline chipSeq pour RNApolII et condensine vont très souvent posséder une portions communes, d'où les bons resultats de l'intersection ! Mais pour autant les positions médianes des pics RNApolII et condensine considérés par HAWKES n'en restent pas moins très éloigné, ce qui explique la vallée en 0 dans l'analyse spatiale.**

**La colocalisation partielle des pics condensine avec ceux de RNApolII peut être interprétée comme la résultante d'un phénomène de compétition entre condensine et RNApolII.**

**On notera que l'utilisation des pics RNApolII Sutani ou Larochelle ne change pas fondamentalement les résultats**

A noter:
*Sutani n'a pas utilisé les données de Kakui pour analyser condensine mais ses propre données*  
La même analyse Hawkes que précédement mais avec condensine Sutani (à la place de condensine Kakui), révèle aussi une vallée en 0. La conclusion de Sutani n'est donc pas liée aux bruits associés aux données SOLID de SUTANI.

![NDRs_union_SRR1564402_rnaPolII_condensine_sutani_merged_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/0faf043121f75b627e4936dad8b916a4/NDRs_union_SRR1564402_rnaPolII_condensine_sutani_merged_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)


## Condensine/NDRs/polII vs TSS/TTS

**N.B: On définit les région TSS/TTS à +/- 150 bases autour du start et stop des gènes connus de S.pombe**  

-> Intersection   

|       Intersect       |          TSS          |     TTS      |
|:---------------------:|:---------------------:|:------------:|
|     Sebastien         |        2726           |     2717     |
|     NDRs union        |        3649           |     3640     |
|     NDRs intersection |      2457             |      2481    |
|     Condensine Kakui  |         248           |     255      |
|    Condensine Sutani  |         187           |     189      |
|condensine intersect   |           102         |    102       |
|condensine union       |        333            |    342       |
| Sutani RNApolII       |         2115          |     2127     |
| Larochelle RNApolII   |         950           |     939      |


**Les NDRs, condensine, et polII semble assez largement colocaliser avec TSS/TTS.**
**A noter une légère diminution de la colocalisation des comptage polII et condensine au TTS par rapport au TSS**

-> analyse spatiale avec RNApol Sutani

**N.B: Seulement l'union des NDRs a été utilisée ici dans l'analyse spatiale**  

rouge: Condensine vs TSS/TTS  
vert: NDRs vs TSS/TTS  
bleu: RNApol vs TSS/TTS

    - TSS
![NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_s_pombe_TSS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/6db9167625af136216dc41bdf99f0356/NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_s_pombe_TSS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

    - TTS

![NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_s_pombe_TTS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg](/uploads/8964e3dd344d31db5a94b03b20ce260a/NDRs_union_SRR1564402_rnaPolII_SRR5227706_condensine_kakui_s_pombe_TTS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_custom.svg)

-> analyse spatiale avec RNApol Sutani

**N.B: Seulement l'union des NDRs a été utilisée ici dans l'analyse spatiale**  

rouge: Condensine vs TSS/TTS  
vert: NDRs vs TSS/TTS  
bleu: RNApol vs TSS/TTS

    - TSS
![NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_s_pombe_TSS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_cutsom.svg](/uploads/b56768fd397d2482f48ee0d4ee5b8ecb/NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_s_pombe_TSS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_cutsom.svg)

    - TTS

![NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_s_pombe_TTS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_cutsom.svg](/uploads/5520d372aa34495dc0ae01bb240839c0/NDRs_union_SRR5227706_condensine_kakui_SRR7289786_larochell_wt_s_pombe_TTS_K_20_delta_1000_kernel_heterogeneous_interval_lambda_1_cutsom.svg)



**Pattern très clair pour les NDRs en amont et en aval du TSS et TTS, respectivement.**  

**Pattern aussi visible pour polII, que ce soit pour Larochelle ou Sutani**
**-> Encore une fois le choix des données n'influence pas fondamentalement les resultats**  

**Plus flou pour condensine**  

## ZOOM rDNA

blanc = gene  
gris/vert = NDRs
rouge = condensine  
bleu= RNApol Larochelle

**N.B: pas RNApolII Sutani dans cette région**

![rDNA.svg](/uploads/0c7a8200773bb8190f576d8e8fee480f/rDNA.svg)

## Résultats analyse de données ChIP-quantitatif

- Mapping

    -> SRR7289786 (IP wt)

    ![reads_origin_SRR7289786.svg](/uploads/3449e3c351db72c27979564b79348833/reads_origin_SRR7289786.svg)

    -> SRR7289787 (WCE wt)

    ![reads_origin_SRR7289787.svg](/uploads/154cea903dfa2d0633d79cf9a9bb8264/reads_origin_SRR7289787.svg)

    -> SRR7289788 (IP mut)

    ![reads_origin_SRR7289788.svg](/uploads/6aef75ae8e44146707d17ae11984d85d/reads_origin_SRR7289788.svg)

    -> SRR7289789 (WCE mut)

    ![reads_origin_SRR7289789.svg](/uploads/99c94dc7e15f5f3fe3b1daf1436c3088/reads_origin_SRR7289789.svg)

- Comparaison IPCwt vs IPCmut  

*N.B: occupancy = couverture par base*

![density_ratio_occupancy.svg](/uploads/cfbbf97f27877d4607e18babe40f45d2/density_ratio_occupancy.svg)

**Les données bruts wt vs mut sont déjà très homogènes !**

- Normalisation : IPX x WCEC / IPC x WCEX

![density_ratio_occupancy_norm.svg](/uploads/a327ac8b5d0b65dff2e61e251f3042cc/density_ratio_occupancy_norm.svg)

**Amélioration du ratio des moyenne d'occupation et diminution de la queue de la distribution à droite.**

- Remove duplicates

![density_ratio_occupancy_dedup_norm.svg](/uploads/56c654f7a907ba383a4536a7ae9df057/density_ratio_occupancy_dedup_norm.svg)

**Le ratio des moyennes d'occupation augmente légèrement mais encore diminution de la queue de la distribution à droite.**
