# Chip analysis pipeline


## Introduction

Nextflow handles job submissions, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

## Reproducibility

In order to guarantee analysis reproducibility, our pipeline use docker images to run bioinfomatics tools.
Note that, each profile configuration used the same tools versions: you will obtained same results, just choose your profile configation depending on container manager available, and the type of machine used.

## Running the pipeline

The typical command for running the pipeline is as follows:
```
./nextflow run src/chip_analysis.nf -c src/chip_analysis.config -profile psmn --inputFromSRA 'data/fromSRA.csv'
--inputFromPATH "data/fromSRA.csv"
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

![CHIP](/uploads/8be6e1185180264cfb5d7ca649653681/CHIP.jpg)

## Main arguments

### -c

As precised above, this pipeline is build to used docker container. The configuration setting needed by our pipeline are precised in **src/chip_analysis.config**

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

Use this to specify the location of csv file . For example:

```
--inputFromPath 'data/fromPATH.csv'
```

Please note the following requirements:

    1) The path of csv file must be enclosed in quotes
    2) CSV file must follow this design

    |           ip            |           input            |
    |:-----------------------:|:--------------------------:|
    | data/sampleIP.fq        | data/sampleINPUT.fq        |      You can specify as many row as you want
    | data/sampleIP_R{1,2}.fq | data/sampleINPUT_R{1,2}.fq |

    3) When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs

    4) Very important: csv field separator must be tabulation

### --inputFromSRA

Use this to specify the SRA ids of your input FastQ files. For example:

```
--inputFromPath 'data/fromSRA.csv'
```

Please note the following requirements:

    1) The path of csv file must be enclosed in quotes
    2) CSV file must follow this design

    |   ip        |  input     |
    |:-----------:|:----------:|          You can specify as many row as you want
    | SRR5227706  | SRR5227705 |

    3) Very important: csv field separator must be tabulation


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

### --macs_gsize

In order to run, MAC2 requires the effective genome size or the size of the genome that is mappable. Mappability is related to the uniqueness of the k-mers at a particular position the genome. Low-complexity and repetitive regions have low uniqueness, which means low mappability. Therefore we need to provide the effective genome length to correct for the loss of true signals in low-mappable regions.
To run macs2 in the pipeline, you have to specify effective length of your genome.

### -resume

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

### --skipFastqc

By default, our pipeline give you the possibility to control the quality of your reads (raw, trimmed our without adapters), by running FastqC. With this option you desactivated FastqC process.

### --skipMultiqc

By default, MultiqC run to summarises all tools reports suitable with multiqc (fastqC, Bowtie{1,2}, cutadapt, picardtools, macs2). This allow you to check user friendly a part of your results.
With this option you desactivated this step.


## Alignment tool

By default, the pipeline uses BOWTIE2 to index fasta genome file, and align the raw FastQ reads to the reference genome.
But if you prefer use BOWTIE1, **--shortReads** option is available.

## MACS2 for peak calling

A commonly used tool for identifying transcription factor binding sites is named Model-based Analysis of ChIP-seq (MACS). The MACS algorithm captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used either for the ChIP sample alone, or along with a control sample which increases specificity of the peak calls.  
Note that our pipeline always needed to provide control/input data linked with IP data. That's why you must follow this design in your CSV input file:

|   ip        |  input     |
|:-----------:|:----------:|     
| SRR5227706  | SRR5227705 |


## Output

Here we describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.
Pipeline overview

The pipeline is built using Nextflow and processes data using the following steps:

    - Cutadapt - adapter trimming
    - Cutadapt - reads trimming
    - FastQC - read quality control
    - Bowtie2/Bowtie1 - alignment
    - Picard tools - remove duplicates
    - MACS2 - peak calling
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

The src/chip_analysis.nf pipeline uses Cutadapt for removal of adapter contamination and trimming of low quality regions. FastQC run after it finishes.

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


### Peak Calling.

MACS2 have been developed for identifying transcription factor binding sites.

```
Output directory: results/peakCalling/sampleIP
                    sampleIP_peak.narrowPeak. #stats, and positions of peaks
                    sampleIP_peaks.xls  #Command line use and summary of peaks
                    sampleIP_macs2_report.txt  #peaks calling stats
```

For further reading and documentation see the MACS2 documentation.


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

For more investigations, follow instructions in *src/hawkes/hawkes* to install and next use *src/hawkes/run_hawkes.R*.

For more informations:
*https://github.com/lereldarion/hawkes*
*doc/biblio/BONNET_2020.pdf*

## Author

* **Arnaud Duvermy**
