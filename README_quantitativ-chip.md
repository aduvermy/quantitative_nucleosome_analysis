# Chip analysis pipeline


## Introduction

Nextflow handles job submissions, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

## Reproducibility

In order to guarantee analysis reproducibility, our pipeline use docker images to run bioinfomatics tools.
Note that, each profile configuration used the same tools versions: you will obtained same results, just choose your profile configation depending on container manager available, and the type of machine used.

## Running the pipeline

The typical command for running the pipeline is as follows:
```
./nextflow run src/chip_quant_analysis.nf -c src/chip_quant_analysis.config -profile psmn --inputFromSRA 'data/fromSRA.csv' --inputFromPATH "data/fromSRA.csv" --macs_gsize "1.25e7"```
This will launch the pipeline with the psmn configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:
```
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
## Workflow

![quantitativ-chip](/uploads/c20b613ad4a3463d4435434443f37392/quantitativ-chip.png)

## Main arguments

### -c

As precised above, this pipeline is build to used docker container. The configuration setting needed by our pipeline are precised in **src/chip_quant_analysis.config**

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

    |           IP_w          |           WCE_w          |           IP_m          |           WCE_m          |
    |:-----------------------:|:--------------------------:|:-----------------------:|:--------------------------:|
    | data/sampleIP_w.fq      | data/sampleINPUT_w.fq      | data/sampleIP_m.fq      | data/sampleINPUT_m.fq      |   
    |data/sampleIP_w_R{1,2}.fq|data/sampleINPUT_w_R{1,2}.fq|data/sampleIP_m_R{1,2}.fq|data/sampleINPUT_m_R{1,2}.fq|

    You can specify as many row as you want.  
    w designe wild type, and m mutant. In that way our pipeline will allow you possibility to compared IP_w with IP_m.  

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

    |           IP_w          |           WCE_w          |           IP_m          |           WCE_m          |
    |:-----------------------:|:--------------------------:|:-----------------------:|:--------------------------:|
    | SRR7289786              | SRR7289787                 | SRR7289788              | SRR7289789                 |

    w designe wild type, and m mutant. In that way our pipeline will allow you possibility to compared IP_w with IP_m.  


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
Follow macs2 

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


## --skipMapping_stats

Our pipeline allow you to visualised part of the reads considered from FASTA or from FASTA CALIBRATION or from both by the mapping step. If you precised **--skipMapping_stats** you will skip this step.
Note that this step is usefull to check quantities of each substances introduced in test tube.

## --skipNorm_stats

Our pipeline allow you to appreciate effect of normalized. If you precised **--skipMapping_stats** you will skip this step.
Note that this step to check if normalized work well on your data.

## MACS2 for peak calling

A commonly used tool for identifying transcription factor binding sites is named Model-based Analysis of ChIP-seq (MACS). The MACS algorithm captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used either for the ChIP sample alone, or along with a control sample which increases specificity of the peak calls.
Note that our pipeline always needed to provide control/input data linked with IP data.


## Output

Here we describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.
Pipeline overview

The pipeline is built using Nextflow and processes data using the following steps:

    - Cutadapt - adapter trimming
    - Cutadapt - reads trimming
    - FastQC - read quality control
    - Bowtie2/Bowtie1 - alignment
    - diagramm_Venn.R - Mappping statistics
    - Picard tools - remove duplicates
    - Deeptools - BAM coverage/occupancy
    - get_normalization_factor.R - Get global normalization factor
    - normalized_occupancy.R - apply normalization factor on each position coverage
    - occupancy_statistics.R - Normalization statistics
    - Samtools - MERGING BAM, Merge IP (wt & mut) and INPUT (wt & mut)
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

The src/chip_quant_analysis.nf pipeline uses Cutadapt for removal of adapter contamination and trimming of low quality regions. FastQC run after it finishes.

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


### Mapping Statistics

We build a R script allow you to directly visualize reads proportion from FASTA or FASTA CALIBRATION or BOTH in a venn diagramm.

```
Output directory: results/mapping/stats
                      sample_reads_origin.svg #Venn diagramm plot

```
Just below, you have an example of output. You can see same quantities of substances from FASTA & from FASTA CALIB have been use in WT and MUT experiment (around 9:1)

![output_mapp_stats](/uploads/3485420053bd89e35f36026e10369a95/output_mapp_stats.png)


### PicardTools

PicardTools locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

```
Output directory: results/mapping/ddup
                      sample_dedup.bam #Compressed alignment file without duplicates
                      sample_dedup_metrics.txt #PicardTools log file
```
For further reading and documentation see the PicardTools documentation.


### Deeptools

This tool takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output. The coverage is calculated as the number of reads per bin, where bins are short consecutive counting windows of a defined size.

```
Output directory: results/bigWig/IP_SAMPLE_id/
                      sample_IP_fasta_occupancy.bw #occupancy norm from fasta
                      sample_IP_fastaCalib_occupancy.bw #occupancy norm from fastaCalib

```

For further reading and documentation see the Deeptools documentation.


### Normalization Statistics

We build a R script allow you to directly visualize effect of normalization on IP calibration.

```
Output directory: results/bigWig/IP_SAMPLE_id/norm
                      sample_IP_fasta_occupancy_norm.bw #occupancy norm from fasta
                      sample_IP_fastaCalib_occupancy_norm.bw #occupancy norm from fastaCalib

```
Just below, an example of normalization effect. On normalized plot, you can see a mean clother to 1, than one you can see on raw plot. More over you can appreciate the decrease of rigth tail distribution. Good news, we improved the fit between IP calib WT & IP calib MUT : normalization work well !

![output_norm_stats](/uploads/51d81ef7a236b1b67a8bca7014ce7f42/output_norm_stats.png)

### Peak Calling

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



## Author

* **Arnaud Duvermy**
