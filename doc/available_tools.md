## Available tools

- **nf module**: a working example of nextflow process is available in `src/nf_modules/<tools>/<tool>.nf` and `src/nf_modules/<tools>/<tool>.config`
- **docker module**: you can create a docker with the `src/docker_modules/<tool>/<version>/docker_init.sh`
- **psmn module**: you can use the tool in the PSMN
- **IN2P3 module**: you can use the tool in the CCIN2P3

| tool | nf module | docker module | psmn module | in2p3 module |
|------|:---------:|:-------------:|:-----------:|:-------------:|
BEDtools | ok | ok | ok | ok 
BFCtools |**no**  | ok | ok | ok
bioawk |**no**  | ok | ok | ok
Bowtie | ok | ok | **no** | ok
Bowtie2 | ok | ok | ok | ok
BWA | ok | ok | ok | ok
canu | ok | ok | ok | ok
cutadapt | ok | ok | ok | ok
deepTools | ok | ok | ok | ok
fastp | ok | ok | ok | ok
FastQC | ok | ok | ok | ok
file_handle | **no** | ok | ok | ok
GATK | **no** | ok | ok | ok
HISAT2 | ok | ok | ok | ok
HTSeq | ok | ok | ok | ok
Kallisto | ok | ok | ok | ok
MACS2 | ok | ok | ok | ok
MultiQC | ok | ok | ok | ok
MUSIC | ok | ok | ok | ok
picard | **no** | ok | ok | ok
pigz | **no** | ok | ok | ok
RSEM | ok | ok | ok | ok
Salmon | **no** | ok | ok | ok
sambamba | ok | ok | ok | ok
samblaster | ok | ok | ok | ok
SAMtools | ok | ok | ok | ok
SRAtoolkit | ok | ok | ok | ok
STAR | ok | ok | ok | ok
subread | **no** | ok | ok | ok
TopHat | **no** | ok | ok | ok
Trimmomatic | **no** | ok | ok | ok
UMItools  | **no** | ok | ok | ok
UrQt | ok | ok | ok | ok
