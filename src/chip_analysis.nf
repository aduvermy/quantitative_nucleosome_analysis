#!/usr/bin/env nextflow

/*
========================================================================================
                         chip_analysis
========================================================================================
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run src/chip_analysis.nf -c src/chip_analysis.config --inputFromSRA data/design/datafromSRA.csv --inputFromPATH data/design/datafromPATH.csv --fasta <genome file> -profile singularity


    Required arguments:
      --inputFromSRA                Full path to directory of CSV file which specifies fastq files SRA IDs(required if --fromPATH not specified)
      --inputFromPATH               Full path to directory of CSV file which specifies full path fastq files: (required if --fromSRA not specified)

    Reference genome
      --fasta                       Full path to directory containing genome fasta file

    Mapping option:
      --shortReads                  Specifies that all input file are not long reads in order to use more adapted mapper: bowtie (default: bowtie2)

    QC Option:
      --skipFastqc                  Skip fastq quality control step (default: activated).
      --skipMultiqc                 Skip merging tools report suitable with multiqc (default: activated)

    Trimming option:
      --trimming                    Activated trimming step (default: desactivated)
      --adapter_removal             Activated adapter removal step (default: desactivated)

    Remove duplicates
      --duplicate_removal           Activated reads duplicates removal step (default: desactivated)

    Peak Calling :
      --macs_gsize                  Effective genome size (if not specified macs2 will not run)

    Nextflow config:
      -c                            Path to config file: src/chip_analysis.config
      -profile                      Profil used by nextflow to run the pipeline (you have choice between singularity, docker, psmn or ccin2p3)
                                    For local utilisation use singularity or docker
    Save option:
      --outdir                      Specify where to save the output from the nextflow run (default: "./results/")

    help message:
      --help                        Print help message
    """
      .stripIndent()
  }


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.help = false
params.trimming = false
params.duplicate_removal = false
params.adapter_removal = false
params.skipFastqc = false
params.skipMultiqc = false
params.shortReads = false
params.inputFromPATH = false
params.inputFromSRA = false
params.fasta = false
params.macs_gsize = false   //for pombe : "1.25e7"
params.outdir = 'results'


/*
 * SET UP CONFIGURATION VARIABLES
 */
//params.help="False"
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Genome']                 = params.fasta ?: 'Not supplied'
summary['INPUT FROM SRA']         = params.inputFromSRA ? params.inputFromSRA : 'Not supplied'
summary['INPUT FROM PATH']        = params.inputFromPATH ? params.inputFromPATH : 'Not supplied'
summary['Trimming Step']          = params.trimming ? 'Yes' :'Skipped'
summary['Remove adapter']         = params.adapter_removal ? 'Yes' : 'Skipped'
summary['Reads QC']               = params.skipFastqc ? 'Skipped' : 'Yes'
summary['Merging Reports']        = params.skipMultiqc ? 'Skipped' : 'Yes'
summary['Mapper']                 = params.shortReads ? 'Bowtie1' : 'Bowtie2'
summary['Remove Duplicate']       = params.duplicate_removal ? 'Yes': 'Skipped'
summary['MACS2 Genome Size']      = params.macs_gsize ?: 'Not supplied'
summary['Config Profile']         = workflow.profile
summary['Output']                 = params.outdir
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       VALIDATE INPUT                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Show a big warning message if we're not running MACS
if (!params.macs_gsize) {
    log.warn "=================================================================\n" +
             "  WARNING! MACS genome size parameter not precised.\n" +
             "  Peak calling analysis will be skipped.\n" +
             "  Please specify value for '--macs_gsize' to run these steps.\n" +
             "======================================================================="
}


if (!params.fasta) {
  exit 1,
    log.warn "=================================================================\n" +
             "  WARNING! No genome fasta file precised.\n" +
             "  Use '--fasta'\n" +
             "  Or '--help' for more informations"
             "======================================================================="
}
else{
    Channel
          .fromPath( params.fasta )
          .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
          .set { fasta_file }
}

if (!params.inputFromSRA && !params.inputFromPATH) {
  exit 1,
   log.warn "=================================================================\n" +
            "  WARNING! No csv input file precised.\n" +
            "  Use '--inputFromSRA' / '--inputFromPATH'" +
            "  Or '--help' for more informations"
            "======================================================================="
}


if (params.inputFromSRA)  {
      Channel
            .fromPath(params.inputFromSRA, checkIfExists: true)
            .set{ fromSRA_csv }
}
else {
      Channel
            .empty()
            .set { fromSRA_csv }
      }

if (params.inputFromPATH)  {
      Channel
            .fromPath(params.inputFromPATH, checkIfExists: true)
            .set { fromPATH_csv }
}
else {
      Channel
            .empty()
            .set { fromPATH_csv }
      }

/*
 *CONCAT INPUT INTO A SINGLE CHANNEL WITH CSV FILE TO CHECKED FORMAT
 * [ fromPATH_csv, fromSRA_csv ]
 */

fromPATH_csv
              .concat(fromSRA_csv)
              .set{ input_csv }



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --            PREPROCESSING - CHECKFORMAT DESIGN INPUT FILE            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process CheckInput_Design {
  tag "$csv.simpleName"
  label "check_input"

  input:
  file csv from input_csv

  output:
  file "*_checked.csv" into design_FromPATH_checked, design_FromSRA_checked
  //python /scratch/Bio/aduvermy/quantitative-nucleosome-analysis/src/check_design.py $design ${design.baseName}_checked.csv
script:
"""
check_design_chip.py $csv ${csv.baseName}_checked.csv
"""
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --            USE CSV INPUT CHECKED TO GET BACK FASTQ FILES            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/*
 * SPLIT DESIGN INPUT CHECKED CHANNEL TO GET BACK FASTQ FILES WITH ADAPTED METHOD
 * [ fromPATH_checked_csv, fromSRA_checked_csv ] -> fromPATH_checked_csv
 *                                               -> fromSRA_checked_csv
 */
if (params.inputFromSRA && params.inputFromPATH){
      design_FromPATH_checked
                              .collect()
                              .first()
                              .set{fromPATH_checked_csv}
      design_FromSRA_checked
                            .collect()
                            .last()
                            .set{fromSRA_checked_csv}
}

if (!params.inputFromSRA){
    design_FromPATH_checked.set{fromPATH_checked_csv}
}

if (!params.inputFromPATH){
    design_FromSRA_checked.set{fromSRA_checked_csv}
}

/*
* CREATE A CHANNEL FOR EXPERIMENT FROM SRA
* [ [ [ file_id_ip, [fastq_files_ip] ], [ file_id_input, [fastq_files_input] ] ] , ... ]
*/

if (params.inputFromSRA){
            fromSRA_checked_csv
                  .splitCsv(header:true, sep:"\t")
                  .map{ row-> tuple(  Channel
                                             .fromSRA(row.ip, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                             .ifEmpty { error "Cannot find any ids matching: ${row.ip}" }
                                             .getVal()
                                        , Channel
                                                .fromSRA(row.input, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                                .ifEmpty { error "Cannot find any ids matching: ${row.input}" }
                                                .getVal()
                                                         ) }
               .set { exp_fromSRA }
}
else{
          Channel
                .empty()
                .set { exp_fromSRA }
}


/*
* CREATE A CHANNEL FOR EXPERIMENT FROM PATH
* [ [ [ file_id_ip, [fastq_files_ip] ], [ file_id_input, [fastq_files_input] ] ] , ... ]
*/

if (params.inputFromPATH){
            fromPATH_checked_csv
                  .splitCsv(header:true, sep:"\t")
                  //.map{row -> tuple (file(row.ip, checkIfExists: true) , file(row.input, checkIfExists: true) ) }
                  .map{ row-> tuple(  Channel
                                             .fromFilePairs(row.ip, size:-1)
                                             .ifEmpty { error "Cannot find any file matching: ${row.ip}" }
                                             .getVal()
                                        , Channel
                                                .fromFilePairs(row.input, size:-1)
                                                .ifEmpty { error "Cannot find any file matching: ${row.input}" }
                                                .getVal()
                                                         ) }
                  .set { exp_fromPATH }
 }
 else{
           Channel
                 .empty()
                 .set { exp_fromPATH }
 }


/*
* CONCAT EXPERIMENT fromPATH AND EXPERIMENT fromSRA CHANNEL INTO A SINGLE EXPERIMENT CHANNEL
* [ [ [ file_id_ip, [fastq_files_ip] ], [ file_id_input, [fastq_files_input] ] ] , ... ]
*/
exp_fromPATH
              .concat(exp_fromSRA)
              .into{ input_ch_2expID; input_ch_2listFQ }


/*
* CREATE CHANNEL FOR EXPERIMENT ID IN ORDER TO MULTIPLEX BAM AFTER
* [ [file_id_ip, file_id_input] , [file_id_ip, file_id_input], ... ]
*/

 input_ch_2expID
       .map { sample -> [ sample[0][0], sample[1][0] ] }
       .set{ experiment_ID }

/*
 * CREATE CHANNEL WITH UNIQ FASTQ FILE
 * [ [ file_id, [fastq_files] ], [ file_id, [fastq_files] ], ... ]
 */

 input_ch_2listFQ
         .flatMap{ sample -> [ sample[0], sample[1] ] }
         .unique()
         .into{ fastq_raw_2QC ; fastq_raw_2trim }


 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 /* --                                                                     -- */
 /* --                         TRIMMING READS                              -- */
 /* --                                                                     -- */
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////

 if (!params.trimming) {
      fastq_raw_2trim.set{ fastq_trim_2cutAdapt }
      Channel
             .empty()
             .set { fastq_trim_2QC }

      Channel
             .empty()
             .set { trimming_report }
   }
 else {
     process Trimming {
       label "cutadapt"
       tag "$file_id"
       publishDir "${params.outdir}/fastq/trim/", mode: 'copy'

       input:
       set file_id, file(reads) from fastq_raw_2trim

       output:
       set file_id, "${file_id}*_trim.fastq.gz" into fastq_trim_2QC, fastq_trim_2cutAdapt
       set file_id, "${file_id}*_report.txt" into trimming_report

       script:
       def single = reads instanceof Path

       if (single){
       """
         cutadapt -q 20,20 \
         -o ${file_id}_trim.fastq.gz \
         ${reads} > ${file_id}_report.txt
       """
       }
       else{
       """
       cutadapt -q 20,20 \
       -o ${file_id}_R1_trim.fastq.gz -p ${file_id}_R2_trim.fastq.gz \
       !{reads[0]} ${reads[1]} > ${file_id}_report.txt
       """
       }
     }
 }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         REMOVE ADAPTERS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

if (!params.adapter_removal) {
      fastq_trim_2cutAdapt.set{ fastq_2mapping}
      Channel
             .empty()
             .set { fastq_cutAdapt_2QC }
      Channel
            .empty()
            .set { adapter_removal_report }
  }
else {
    process Adapter_removal {
      label "cutadapt"
      tag "$file_id"
      publishDir "${params.outdir}/fastq/cut/", mode: 'copy'

      input:
      set file_id, file(reads) from fastq_trim_2cutAdapt

      output:
      set file_id, "${file_id}*_cut.fastq.gz" into fastq_2mapping, fastq_cutAdapt_2QC
      set file_id, "${file_id}_report.txt" into adapter_removal_report

      script:
      def single = reads instanceof Path

      if (single){
      """
      cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT\
      -o ${file_id}_cut.fastq.gz \
      ${reads} > ${file_id}_report.txt
      """
      }
      else{
      """
      cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
      -o ${file_id}_R1_cut.fastq.gz -p ${file_id}_R2_cut.fastq.gz \
      ${reads[0]} ${reads[1]} > ${file_id}_report.txt
      """
      }
}
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       READS QUALITY CONTROLE                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/*
 * CONCAT RAW, TRIM, CUTADAPT FASTQ INTO SINGLE Channel
 */
fastq_raw_2QC
            .concat(fastq_trim_2QC)
            .concat(fastq_cutAdapt_2QC)
            .set{fastq_2QC}


if (params.skipFastqc) {
        Channel
              .empty()
              .set { fastqc_report }
}
else{
process Fastqc {
     label "fastqc"
     tag "$file_id"
     publishDir "${params.outdir}/fastq/QC/", mode: 'copy'

     input:
     set file_id, file(reads) from fastq_2QC

     output:
     file "*.{zip,html}" into fastqc_report

     script:
     def single = reads instanceof Path

     if (single){
        """
        fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
        ${reads}
        """
      }
    else{
        """
        fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
        ${reads[0]} ${reads[1]}
        """
      }
}
}






///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   UNZIP FASTA FILE IF NEEDED                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * GET EXTENSION FASTA FILE TO TEST IT
 */

fasta_file.into{ fasta_file_test_zip;
                 fasta_file_zip }

ext=fasta_file_test_zip.getVal().getExtension()


 if(ext=="gz" || ext=="bz" || ext=="zip"){

    log.info "Genome fasta file zip, we need to unzip it for next steps"

    process Unzip_fasta {
     tag "$fasta.simpleName"

     input:
       file fasta from fasta_file_zip

     output:
       file "${fasta.simpleName}.fasta" into fasta_2indexing

     script:
     """
     zcat ${fasta} > ${fasta.simpleName}.fasta
     """
    }
}
else {
  fasta_file_zip.set { fasta_2indexing }
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      BUILD GENOME INDEX                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process Fasta_indexing {
    if (!params.shortReads){
          label "bowtie2"
    }
    else {
          label "bowtie"
    }
    tag "$fasta.simpleName"

    input:
    file fasta from fasta_2indexing

    output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

    script:
    if (params.shortReads){
    """
    bowtie-build --threads ${task.cpus} -f ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie_report.txt

    if grep -q "Error" ${fasta.baseName}_bowtie_report.txt; then
      exit 1
    fi
    """
    }
    else{
    """
    bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie2_report.txt

    if grep -q "Error" ${fasta.baseName}_bowtie2_report.txt; then
      exit 1
    fi
    """
    }
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MAPPING                                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process Mapping {
   publishDir "${params.outdir}/mapping/${file_id}", mode: 'copy'

   if (!params.shortReads){
     label "bowtie2"
   }
   else {
     label "bowtie"

   }

  tag "$file_id"

   input:
   set file_id, file(reads) from fastq_2mapping
   file index from index_files.collect()


   output:
   set file_id, val(data_type) ,"*.bam" into bam_files
   file "*_report.txt" into mapping_report

   script:
   def single = reads instanceof Path

   if (params.shortReads){
     index_id = index[0]
     for (index_file in index) {
       if (index_file =~ /.*\.1\.ebwt/ && !(index_file =~ /.*\.rev\.1\.ebwt/)) {
           index_id = ( index_file =~ /(.*)\.1\.ebwt/)[0][1]
         }
     }
     if (single){
       data_type="SE"
       """
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       -q ${reads} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report.txt
       """
     }
     else{
       data_type="PE"
       """
       # -v specify the max number of missmatch, -k the number of match reported per
       # reads
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}.bam


       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report.txt
       """

   }
 }
   else {
     index_id = index[0]
     for (index_file in index) {
       if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
           index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
       }
     }
     if (single){
       data_type="SE"
       """
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       -U ${reads} 2> \
       ${file_id}_bowtie2_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}.bam

       if grep -q "Error" ${file_id}_bowtie2_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie2_report_tmp.txt > ${file_id}_bowtie2_report.txt
       """
     }
     else{
       data_type="PE"
       """
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie2_report_tmp.txt | \
       samtools view -Sbh -@ ${task.cpus} - | \
       samtools sort -@ ${task.cpus} -o ${file_id}.bam

       if grep -q "Error" ${file_id}_bowtie2_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie2_report_tmp.txt > ${file_id}_bowtie2_report.txt
       """

   }
}
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       REMOVE DUPLICATES                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

if (!params.duplicate_removal){
  dedup_bam_files=bam_files
  Channel
        .empty()
        .set { dDup_bam_report }

}
else{
    process Duplicate_removal {
      label 'picardtools'
      tag "$file_id"
      publishDir "${params.outdir}/mapping/${file_id}/ddup/", mode: 'copy'

      input:
        set file_id, data_type, file(bam) from bam_files
      output:
        set file_id , data_type ,"*_dedup.bam" into dedup_bam_files
        set file_id, "*.txt" into dDup_bam_report

      script:
      """
      PicardCommandLine MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${bam} \
      OUTPUT=${file_id}_dedup.bam \
      METRICS_FILE=${file_id}_metrics.txt
      """
      }
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         PEAK CALLING                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/*
 * MULTIPLEXAGE, AFFECT BAM TO EXPERIMENT(IP,INPUT)
 * [ [ip_ID, input_ID, bam_ip, bam_input, data_type ], ...]
 */

dedup_bam_files.into{ bam2match_ip ; bam2match_input}

experiment_ID
  .combine( bam2match_ip, by:0 )
  .combine( bam2match_input.map{ it-> [ it[1], it[0],it[2] ] }, by:1)
  .map{ it  -> [ it[1] , it[0] , it[3] , it[5], it[2] ] }
  .set{ ip_vs_input }


process Peak_calling {
  tag "${ip_ID} vs ${input_ID}"
  label "macs2"
  publishDir "${params.outdir}/peakCalling/${ip_ID}", mode: 'copy'

  input:
    set ip_ID, input_ID, file(bam_ip), file(bam_input), data_type from ip_vs_input

  output:
    file "*" into peak_output
    file "*_peaks.xls" into peak_calling_report


  when:
  params.macs_gsize

  script:
  format = data_type=="SE" ? "BAM" : "BAMPE"
  """
  macs2 callpeak \
    --fix-bimodal \
    -f ${format} \
    --treatment ${bam_ip} \
    --control ${bam_input} \
    --name ${ip_ID} \
    --gsize ${params.macs_gsize} 2> \
    ${ip_ID}_macs2_report.txt

  if grep -q "ERROR" ${ip_ID}_macs2_report.txt; then
    echo "MACS2 error"
    exit 1
  fi
  """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      MERGE ALL STEPS REPORTS                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  process MultiQC {
     label "multiQC"
     publishDir "${params.outdir}/multiQC", mode: 'copy'

     input:
     file report_fastqc from fastqc_report.collect().ifEmpty([])
     file report_trim from trimming_report.collect().ifEmpty([])
     file report_adptoRemoval from adapter_removal_report.collect().ifEmpty([])
     file report_mapping from mapping_report.collect()
     file report_ddup from dDup_bam_report.collect().ifEmpty([])
     file report_peakC from peak_calling_report.collect().ifEmpty([])

     output:
     file "*multiqc_*" into multiqc_report


     when:
     !params.skipMultiqc

     script:
     """
     multiqc -f . \\
     -m fastqc -m cutadapt -m bowtie1 -m bowtie2 -m picard -m macs2
     """
  }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       NF-CORE HEADER                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset =  "\033[0m";
    c_dim = "\033[2m";
    c_black = "\033[0;30m";
    c_green = "\033[0;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_purple = "\033[0;35m";
    c_cyan = "\033[0;36m";
    c_white ="\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  CHIP-SEQ Pipeline          ${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
