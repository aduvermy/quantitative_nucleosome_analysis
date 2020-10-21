# MNASE analysis pipeline


## Introduction

Nextflow handles job submissions, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

## Reproducibility

In order to guarantee analysis reproducibility, our pipeline use docker images to run bioinfomatics tools.
Note that, each profile configuration used the same tools versions: you will obtained same results, just choose your profile configation depending on container manager available, and the type of machine used.

## Running the pipeline

The typical command for running the pipeline is as follows:
```
./nextflow run src/mnase_analysis.nf -c src/mnase_analysis.config -profile psmn --inputFromSRA 'ERR1353027'
--inputFromPATH "data/sample_{1,2}.fq"
```
This will launch the pipeline with the psmn configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:
```
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
## Workflow

![Mnase_pipeline](/uploads/7add10ecb648743eb4caa9981012efd8/Mnase_pipeline.jpg)

## Main arguments

### -c

As precised above, this pipeline is build to used docker container. The configuration setting needed by our pipeline are precised in **src/mnase_analysis.config**

If -c is not specified at all the pipeline will be run locally and expects all software to be installed and available on the PATH, with adapted version.

### -profile

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: -profile psmn


  - psmn  
Use this configuration if you run the pipeline on lbmc cluster. Settings are configured to used resources available.

  - docker  
Local utilisation, if docker is installed on your machine.

  - singularity  
Local utilisation, if singularity is installed on your machine.


### --inputFromPATH

Use this to specify the location of your input FastQ files. For example:

```
--inputFromPath 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

    1) The path must be enclosed in quotes
    2) When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.
    3) You can specify many samples using list


### --inputFromSRA

Use this to specify the SRA ids of your input FastQ files. For example:

```
--inputFromSRA ["ERR1353027", "ERR1353028", "ERR1353029"]
```

Please note the following requirements:

    1) The path must be enclosed in quotes
    3) You can specify many samples using list

**N.B:** *--inputFromPATH* **and** *--inputFromSRA* **can be used** *simultaneously.*  
*Paired end* **and** *single end* **data can be specify simultaneously with both, inputFromSRA and inputFromPATH. Our pipeline will detect number of file involved with each sample (1 or 2), and will adapt treatment for each type of data.**  
**So careful to well follow syntax needed to precised path location of paired end data, as see above.**


### --fasta  

Full path to directory containing genome fasta file, used as reference for mapping.


### --trimming:

By default, raw FastQ reads are used for mapping but you can choose to trimmed your read, with this option.


### --adapter_removal             

You can also choose to remove Illumina adapter in your data, with this option (default: desactivated).

### --duplicate_removal

Reads duplicates can be a source of problems, in bioinfomatics analysis. This option allow you to remove reads duplicates (default: desactivated).


### -resume

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

### --skipFastqc

By default, our pipeline give you the possibility to control the quality of your reads (raw, trimmed our without adapters), by running FastqC. With this option you desactivated FastqC process.

### --skipMultiqc

By default, MultiqC run to summarises all tools reports suitable with multiqc (fastqC, Bowtie{1,2}, cutadapt, picardtools). This allow you to check user friendly a part of your results.  
With this option you desactivated this step.


## Alignment tool

By default, the pipeline uses BOWTIE2 to index fasta genome file, and align the raw FastQ reads to the reference genome.
But if you prefer use BOWTIE1, **--shortReads** option is available.

## DANPOS for peak calling

In this pipeline, peak calling is realized by algorithm **dpos** of the **DANPOS** tools. This tools is well adapted to analysed MNASE-seq data.

For more informations about DANPOS:
*https://sites.google.com/site/danposdoc/*  
*doc/biblio/DANPOS.pdf*


## Peak Calling normalization

Peak calling process will results in wig format file, which count number of reads covering each base of the reference (occupancy).  
In order to locate *nucleosome depleted regions*, occupancy needs to be normalized to prevent amplification bias.  
For this  you can choose to normalized each position by average occupancy along genome (default behaviour), or you can used *quantile normalization* suggest by **DANPOS**, with **--quantile_norm** option.
Note that we not yet got quantile normalization right.

## NDRs Finder

We create a python script to parsed .wig file normalized, and build .bed file identifying regions seems to be NDRs.
For this, we considered as **nucleosome-depleted regions (NDRs) any chromosomal region of at least 150 bp in length exhibiting a normalized MNase-seq coverage depth smaller than 0.4 (Soriano et al, 2013).**  
Note that this definition is well adapted with average occupancy normalization.

## Output

Here we describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.
Pipeline overview

The pipeline is built using Nextflow and processes data using the following steps:

    - Cutadapt - adapter trimming
    - Cutadapt - reads trimming
    - FastQC - read quality control
    - Bowtie2/Bowtie1 - alignment
    - Picard tools - remove duplicates
    - DANPOS - peak calling
    - DANPOS/wig_average_normalization.py - quantile wig normalization/average wig normalization
    - DANPOS - smooth wig normalized
    - wig2bed.py - Find NDRs regions
    - MultiQC - aggregate report, describing results of the whole pipeline

### FastQC

FastQC gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

```
Output directory: results/fastq/QC
                      sample_fastqc.html #FastQC report, containing quality metrics for your fastq files
                      zips/sample_fastqc.zip #zip file containing the FastQC report,
                                             #tab-delimited data file and plot images
```

For further reading and documentation see the FastQC documentation.


### Cutadapt

The src/mnase_analysis.nf pipeline uses Cutadapt for removal of adapter contamination and trimming of low quality regions. FastQC run after it finishes.

MultiQC reports the percentage of bases removed by Cutadapt in the General Statistics table, along with a line plot showing where reads were trimmed.

```
Trimming output directory: results/fastq/trim
Remove adapters output directory: results/fastq/cut
#Contains FastQ files with quality and adapter trimmed reads for each sample, along with a log file describing the trimming.
```
For further reading and documentation see the Cutadapt documentation.

### Bowtie2/Bowtie1

Bowtie is a read aligner designed for data sequencing.

The Bowtie section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as Uniquely mapped and few Unmapped reads.

```
Output directory: results/mapping
                      sample.bam  #Compressed alignment file
                      sample_report.txt #Bowtie1/2 log file
```

For further reading and documentation see the Bowtie documentation.

### PicardTools

PicardTools locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

```
Output directory: results/mapping/ddup
                      sample_dedup.bam #Compressed alignment file without duplicates
                      sample_dedup_metrics.txt #PicardTools log file
```

For further reading and documentation see the PicardTools documentation.

### Peak Calling to locate NDRs

DANPOS have been developed for peak Calling. Dpos algorithm is particulary adapted for analyse nucleosomes positions.

```
Output directory: results/peakCalling/sample
                    sample_positions.xls #stats, and positions of peaks
                    sample.wig  #Contain protein occupancy values at each base pair across the whole genome
```

The wig file is next normalized, depending on the method choose (quantile or average normalization), and then smooth.

```
Output directory: results/peakCalling/sample/normalization
                        sample.{a,q}nor.wig #Contain protein occupancy values normalized at each base pair
                                            #across the whole genome
                        sample.{a,q}nor.wig #Contain protein occupancy values normalized and smooth at each base
                                            #pair across the whole genome
                        sample.{a,q}nor.smooth.positions.xls #stats, and positions of peaks norm and smooth
```

By using **sample.{a,q}nor.wig** file we can now, defined nucleosome depleted region, as defined in Soriano et al (2013).  
```
Output directory: results/peakCalling/sample/NDRs
                      sample.bed  #Regions identify as corresponding to NDRs definition
```
For further reading and documentation see the DANPOS documentation.


## MultiQC

MultiQC is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

```
Output directory: results/multiqc
                        multiqc_report.html #MultiQC report - HTML file that can be viewed in your web browser
                        Project_multiqc_data/ #Directory containing parsed statistics from the
                                              #different tools used in the pipeline
```

For further reading and documentation see the MultiQC documentation.

## Hawkes

for more investigations, follow instructions in *src/hawkes/hawkes* to install and next use *src/hawkes/run_hawkes.R*.

For more informations:
*https://github.com/lereldarion/hawkes*  
*doc/biblio/BONNET_2020.pdf*

## Author

* **Arnaud Duvermy**
