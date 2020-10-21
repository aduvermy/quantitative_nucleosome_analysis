#!/usr/bin/env nextflow

/*
========================================================================================
                         QUANTITATIVE chip_analysis
========================================================================================
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run src/chip_quant_analysis.nf -c src/chip_quant_analysis.config --inputFromSRA datafromSRA.csv --inputFromPATH datafromPATH.csv --fasta <genome file> -profile singularity


    Required arguments:
      --inputFromSRA                Full path to directory of CSV file which specifies fastq files SRA IDs(required if --fromPATH not specified)
      --inputFromPATH               Full path to directory of CSV file which specifies full path fastq files: (required if --fromSRA not specified)

    Reference genome
      --fasta                       Full path to directory containing genome fasta file
      --fasta_calib                 Full path to directory containing genome fasta file used for calibration


    Mapping option:
      --shortReads                  Specifies that all input file are not long reads in order to use more adapted mapper: bowtie (default: bowtie2)

    QC Option:
      --skipFastqc                  Skip fastq quality control step (default: activated).
      --skipMultiqc                 Skip merging tools report suitable with multiqc (default: activated)

    Statistics Option:
      --skipMapping_stats           Skip building Venn diagramm to check reads origin proportion (default: activated)
      --skipNorm_stats              Skip building occupancy plot to check normalization behaviour (default: activated)

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






////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.help = false
params.adapter_removal = false
params.trimming = false
params.skipFastqc = false
params.skipMultiqc = false
params.shortReads = false
params.duplicate_removal = false
params.skipNorm_stats = false
params.skipMapping_stats = false

params.inputFromPATH = false
params.inputFromSRA= false
params.fasta= false
params.fasta_calib= false
params.macs_gsize = false
params.outdir = 'results'



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
summary['Genome calibration']     = params.fasta_calib ?: 'Not supplied'
summary['INPUT FROM SRA']         = params.inputFromSRA ? params.inputFromSRA : 'Not supplied'
summary['INPUT FROM PATH']        = params.inputFromPATH ? params.inputFromPATH : 'Not supplied'
summary['Trimming Step']          = params.trimming ? 'Yes' :'Skipped'
summary['Remove adapter']         = params.adapter_removal ? 'Yes' : 'Skipped'
summary['Reads QC']               = params.skipFastqc ? 'Skipped' : 'Yes'
summary['Merging Reports']        = params.skipMultiqc ? 'Skipped' : 'Yes'
summary['Mapper']                 = params.shortReads ? 'Bowtie1' : 'Bowtie2'
summary['Mapping statistics']     = params.skipMapping_stats ? 'Skipped' : 'Yes'
summary['Remove Duplicate']       = params.duplicate_removal ? 'Yes': 'Skipped'
summary['MACS2 Genome Size']      = params.macs_gsize ?: 'Not supplied'
summary['Normalization stats']    = params.skipNorm_stats ? 'Skipped' : 'Yes'
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
             "  Use --fasta <my genome> "
             "  Or --help for more informations"
             "======================================================================="
}
else{
  Channel
        .fromFilePairs( params.fasta , size:-1)
        .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
        .into { fasta_2concat ; fasta_2joinIndex }
}



if (!params.fasta_calib) {
    exit 1,
    log.warn "=================================================================\n" +
             "  WARNING! No Calibration genome fasta file precised.\n" +
             "  Use '--fasta_calib' \n"
             "  Or '--help' for more informations"
             "======================================================================="
}
else{
     Channel
           .fromFilePairs( params.fasta_calib , size:-1)
           .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta_calib}" }
           .into { fastaCalib_2concat; fastaCalib_2joinIndex }
}




if (!params.inputFromSRA && !params.inputFromPATH) {
   exit 1,
   log.warn "=================================================================\n" +
            "  WARNING! No csv input file precised.\n" +
            "  Use '--inputFromSRA' / '--inputFromPATH' \n" +
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
check_design_chip_quant.py $csv ${csv.baseName}_checked.csv
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
                                             .fromSRA(row.IP_w, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                             .ifEmpty { error "Cannot find any ids matching: ${row.IP_w}" }
                                             .getVal()
                                        , Channel
                                                .fromSRA(row.WCE_w, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                                .ifEmpty { error "Cannot find any ids matching: ${row.WCE_w}" }
                                                .getVal()
                                        ,  Channel
                                                 .fromSRA(row.IP_m, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                                 .ifEmpty { error "Cannot find any ids matching: ${row.IP_m}" }
                                                 .getVal()
                                            , Channel
                                                    .fromSRA(row.WCE_m, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
                                                    .ifEmpty { error "Cannot find any ids matching: ${row.WCE_m}" }
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
                                             .fromFilePairs(row.IP_w, size:-1)
                                             .ifEmpty { error "Cannot find any file matching: ${row.IP_w}" }
                                             .getVal()
                                        , Channel
                                                .fromFilePairs(row.WCE_w, size:-1)
                                                .ifEmpty { error "Cannot find any file matching: ${row.WCE_w}" }
                                                .getVal()
                                        , Channel
                                                 .fromFilePairs(row.IP_m, size:-1)
                                                 .ifEmpty { error "Cannot find any file matching: ${row.IP_m}" }
                                                 .getVal()
                                        , Channel
                                                .fromFilePairs(row.WCE_m, size:-1)
                                                .ifEmpty { error "Cannot find any file matching: ${row.WCE_m}" }
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
* [ [ [ file_id_IP_w, [fastq_files_IP_w] ], [ file_id_WCE_w, [fastq_files_WCE_w] ] ] , ... ]
*/
exp_fromPATH
              .concat(exp_fromSRA)
              .into{ input_ch_2expID; input_ch_2listFQ }


/*
* CREATE CHANNEL FOR EXPERIMENT ID IN ORDER TO MULTIPLEX BAM AFTER
* [ [file_id_IP_w, file_id_WCE_w, file_id_IP_m, file_id_WCE_m] , [file_id_IP_w, file_id_WCE_w, file_id_IP_m, file_id_WCE_m], ... ]
*/

 input_ch_2expID
       .map { sample -> [ sample[0][0], sample[1][0], sample[2][0], sample[3][0] ] }
       .into{ experiment_ID_2bam_WT_MUT ; experiment_ID_2log_2IP_vs_INPUT; experiment_ID_2bw_IP_WT_MUT ; experiment_ID_2bw_norm_IP_WT_MUT }

/*
 * CREATE CHANNEL WITH UNIQ FASTQ FILE
 * [ [ file_id, [fastq_files] ], [ file_id, [fastq_files] ], ... ]
 */

 input_ch_2listFQ
         .flatMap{ sample -> [ sample[0], sample[1], sample[2], sample[3] ] }
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

fasta_2concat
      .concat(fastaCalib_2concat)
      .map{ it -> [ it[0], it[1][0].getExtension(), it[1][0] ] }
      .set{ fasta_2unzip }


process CheckFasta_format {
 tag "$fasta.simpleName"

 input:
   set fasta_id, val(extension), file(fasta)  from fasta_2unzip

 output:
   set fasta_id, "${fasta.simpleName}_unzip.fasta" into fasta_2indexing

 script:
    if (extension=="gz" || extension=="bz" || extension=="zip"){
      log.info "${fasta.simpleName} file zip, we need to unzip it for next steps"
      """
      zcat ${fasta} > ${fasta.simpleName}_unzip.fasta
      """
    }
    else{
      """
      mv ${fasta} ${fasta.simpleName}_unzip.fasta
      """
    }

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
   tag "$fasta_id"

   input:
   set fasta_id , file (fasta) from fasta_2indexing

   output:
   set fasta_id, "*.index*" into index_2joinfasta, index_2joinfastaCalib
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


/*
 *SPLIT FASTA INDEX CH FOR MAPPING
 * [ fasta_calib_file, fasta_file ]
 */


index_2joinfasta
                      .join(fasta_2joinIndex)
                      .map{ it -> it[1]}
                      .set{ index_fasta_files }
index_2joinfastaCalib
                      .join(fastaCalib_2joinIndex)
                      .map{ it -> it[1]}
                      .set{ index_calib_files }




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MAPPING                                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



process Mapping {

   publishDir "${params.outdir}/mapping/${file_id}", mode: 'copy', pattern:"*.bam"
   publishDir "${params.outdir}/fastq/${file_id}", mode: 'copy', pattern:"*.fastq"

   tag "$file_id"

   if (!params.shortReads){
     label "bowtie2"
   }
   else {
     label "bowtie"
   }

   input:
   set file_id, file(reads) from fastq_2mapping
   file index_calib from index_calib_files.collect()
   file index from index_fasta_files.collect()


   output:
   set file_id, val(data_type) ,"*.bam" into bam_files
   set file_id, "*.fastq" into reads
   set file_id, "*_bowtie_report_exclusiv2fasta.txt", "*_bowtie_report_both.txt", "*_bowtie_report_exclusiv2fastaCalib.txt" into mapping_reports2stats, mapping_reports2normalization

   script:
   def single = reads instanceof Path

   if (params.shortReads){
     index_id = index[0]
     for (index_file in index) {
       if (index_file =~ /.*\.1\.ebwt/ && !(index_file =~ /.*\.rev\.1\.ebwt/)) {
           index_id = ( index_file =~ /(.*)\.1\.ebwt/)[0][1]
         }
     }
     index_c_id = index_calib[0]
     for (index_file in index_calib) {
       if (index_file =~ /.*\.1\.ebwt/ && !(index_file =~ /.*\.rev\.1\.ebwt/)) {
           index_c_id = ( index_file =~ /(.*)\.1\.ebwt/)[0][1]
         }
     }
     if (single){
       data_type="SE"
       """
       ##ALIGNED RAW TO FASTA
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       --un ${file_id}_unaligned_2fasta.fastq --al ${file_id}_aligned_2fasta.fastq \
       -q ${reads} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report.txt

       ##ALIGNED RAW TO FASTA CALIB
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --un ${file_id}_unaligned_2fastaCalib.fastq --al ${file_id}_aligned_2fastaCalib.fastq \
       -q ${reads} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_calib.txt

       ##ALIGNED READS, ALREADY MAPPED TO FASTA, TO FASTA CALIB NOW -> get multiMap
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --al ${file_id}_aligned_2Both.fastq \
       -q ${file_id}_aligned_2fasta.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2both.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_both.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA CALIB, TO FASTA  NOW -> get fasta exclusif map
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       --al ${file_id}_exclusiv_aligned_2fasta.fastq \
       -q ${file_id}_unaligned_2fastaCalib.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fasta.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA, TO FASTA CALIB  NOW -> get fasta calib exclusif map
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --al ${file_id}_exclusiv_aligned_2fasta.fastq \
       -q ${file_id}_unaligned_2fasta.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fastaCalib.txt

       ##rm useless files
       rm ${file_id}_2fasta.bam ${file_id}_2fastaCalib.bam ${file_id}_unaligned_2fasta*.fastq ${file_id}_aligned_2fasta*.fastq
       """
     }
     else{
       data_type="PE"
       """
       ##ALIGNED RAW TO FASTA
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       --un ${file_id}_unaligned_2fasta.fastq --al ${file_id}_aligned_2fasta.fastq \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report.txt

       ##ALIGNED RAW TO FASTA CALIB
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --un ${file_id}_unaligned_2fastaCalib.fastq --al ${file_id}_aligned_2fastaCalib.fastq \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_calib.txt

       ##ALIGNED READS, ALREADY MAPPED TO FASTA, TO FASTA CALIB NOW -> get multiMap
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --al ${file_id}_aligned_2Both.fastq \
       -1 ${file_id}_aligned_2fasta_1.fastq -2 ${file_id}_aligned_2fasta_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2both.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_both.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA CALIB, TO FASTA  NOW -> get fasta exclusif map
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
       --al ${file_id}_exclusiv_aligned_2fasta.fastq \
       -1 ${file_id}_unaligned_2fastaCalib_1.fastq -2 ${file_id}_unaligned_2fastaCalib_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fasta.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA, TO FASTA CALIB  NOW -> get fasta calib exclusif map
       bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_c_id} \
       --al ${file_id}_exclusiv_aligned_2fasta.fastq \
       -1 ${file_id}_unaligned_2fasta_1.fastq -2 ${file_id}_unaligned_2fasta_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fastaCalib.txt

       ##rm useless files
       rm ${file_id}_2fasta.bam ${file_id}_2fastaCalib.bam ${file_id}_unaligned_2fasta*.fastq ${file_id}_aligned_2fasta*.fastq
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

       index_c_id = index_calib[0]
       for (index_file in index_calib) {
         if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
             index_c_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
         }
     }
     if (single){
       data_type="SE"
       """
       ##ALIGNED RAW TO FASTA
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       --un ${file_id}_unaligned_2fasta.fastq --al ${file_id}_aligned_2fasta.fastq \
       -U ${reads} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report.txt

       ##ALIGNED RAW TO FASTA CALIB
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --un ${file_id}_unaligned_2fastaCalib.fastq --al ${file_id}_aligned_2fastaCalib.fastq \
       -U ${reads} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_calib.txt

       ##ALIGNED READS, ALREADY MAPPED TO FASTA, TO FASTA CALIB NOW -> get multiMap
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --al ${file_id}_aligned_2Both.fastq \
       -U ${file_id}_aligned_2fasta.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2both.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_both.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA CALIB, TO FASTA  NOW -> get fasta exclusif map
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       --al ${file_id}_exclusiv_aligned_2fasta.fastq \
       -U ${file_id}_unaligned_2fastaCalib.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fasta.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA, TO FASTA CALIB  NOW -> get fasta calib exclusif map
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --al ${file_id}_exclusiv_aligned_2fastaCalib.fastq \
       -U ${file_id}_unaligned_2fasta.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fastaCalib.txt

       ##rm useless files
       rm ${file_id}_2fasta.bam ${file_id}_2fastaCalib.bam ${file_id}_unaligned_2fasta*.fastq ${file_id}_aligned_2fasta*.fastq
       """
     }
     else{
       data_type="PE"
       """
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       --un-conc ${file_id}_unaligned_2fasta_%.fastq --al-conc ${file_id}_aligned_2fasta_%.fastq \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie2_report_tmp.txt | \
       samtools view -Sbh -@ ${task.cpus} - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fasta.bam

       if grep -q "Error" ${file_id}_bowtie2_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie2_report_tmp.txt > ${file_id}_bowtie2_raw2fasta.txt


       ##ALIGNED RAW TO FASTA CALIB
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --un-conc ${file_id}_unaligned_2fastaCalib_%.fastq --al-conc ${file_id}_aligned_2fastaCalib_%.fastq \
       -1 ${reads[0]} -2 ${reads[1]} 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_raw2calib.txt

       ##ALIGNED READS, ALREADY MAPPED TO FASTA, TO FASTA CALIB NOW -> get multiMap
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --al-conc ${file_id}_aligned_2Both_%.fastq \
       -1 ${file_id}_aligned_2fasta_1.fastq -2 ${file_id}_aligned_2fasta_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_2both.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_both.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA CALIB, TO FASTA  NOW -> get fasta exclusif map
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
       --al-conc ${file_id}_exclusiv_aligned_2fasta_%.fastq \
       -1 ${file_id}_unaligned_2fastaCalib_1.fastq -2 ${file_id}_unaligned_2fastaCalib_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fasta.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fasta.txt

       ##ALIGNED READS, ALREADY UNMAPPED TO FASTA, TO FASTA CALIB  NOW -> get fasta calib exclusif map
       bowtie2 --very-sensitive -p ${task.cpus} -x ${index_c_id} \
       --al-conc ${file_id}_exclusiv_aligned_2fastaCalib_%.fastq \
       -1 ${file_id}_unaligned_2fasta_1.fastq -2 ${file_id}_unaligned_2fasta_2.fastq 2> \
       ${file_id}_bowtie_report_tmp.txt | \
       samtools view -@ ${task.cpus} -Sbh - | \
       samtools sort -@ ${task.cpus} -o ${file_id}_exclusiv2fastaCalib.bam

       if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
       exit 1
       fi
       tail -n 19 ${file_id}_bowtie_report_tmp.txt > ${file_id}_bowtie_report_exclusiv2fastaCalib.txt

       ##rm useless files
       rm ${file_id}_2fasta.bam ${file_id}_2fastaCalib.bam ${file_id}_unaligned_2fasta*.fastq ${file_id}_aligned_2fasta*.fastq
       """

   }
}
}





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       MAPPING STATISTICS                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



if (params.skipMapping_stats){
  Channel
        .empty()
        .set { plots_mapping }

}
else{
    process Statistics_mapping {
      label 'R'
     tag "$file_id"
     publishDir "${params.outdir}/mapping/${file_id}/stats", mode: 'copy'

     input:
        set file_id, report2fasta, report2both, report2fastaCalib from mapping_reports2stats

     output:
       file "${file_id}_reads_origin.svg" into plots_mapping

     script:
     if (params.shortReads){
       tool=1
     }
     else{
       tool=2
     }
     """
     diagramm_Venn.R \
     ${report2fasta} \
     ${report2both} \
     ${report2fastaCalib} \
     ${tool} \
     ${file_id}_reads_origin.svg
     """
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
        set file_id, data_type, file(bams) from bam_files
      output:
        set file_id , data_type , "*.bam" into dedup_bam_files
        set file_id, "${bams[1].baseName}_metrics.txt", "${bams[0].baseName}_metrics.txt", "${bams[2].baseName}_metrics.txt" into dDup_bam_report, dedup_reports2normalization


      script:
      """
      PicardCommandLine MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=true \
      INPUT=${bams[0]} \
      OUTPUT=${bams[0].baseName}_dedup.bam \
      METRICS_FILE=${bams[0].baseName}_metrics.txt

      PicardCommandLine MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=true \
      INPUT=${bams[1]} \
      OUTPUT=${bams[1].baseName}_dedup.bam \
      METRICS_FILE=${bams[1].baseName}_metrics.txt

      PicardCommandLine MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=true \
      INPUT=${bams[2]} \
      OUTPUT=${bams[2].baseName}_dedup.bam \
      METRICS_FILE=${bams[2].baseName}_metrics.txt
      """
      }
}



/*
 * MULTIPLEXAGE, AFFECT BAM TO EXPERIMENT(IP,INPUT)
 * [ [ip_ID, input_ID, bam_ip, bam_input, data_type ], ...]
 */

dedup_bam_files.into{ bam2match_ip_wt ; bam2match_input_wt; bam2match_ip_mut ; bam2match_input_mut}


experiment_ID_2bam_WT_MUT
      .combine( bam2match_ip_wt.map{ it-> [ it[0], it[1],[ it[2][1], it[2][2] ] ] }, by:0 )  //match bam with IP_wt
      .combine( bam2match_input_wt.map{ it-> [ it[1], it[0],[ it[2][1], it[2][2] ] ] }, by:1 ) //match bam with INPUT_wt
      .map{ it  -> [ it[2], it[3], [ it[1], it[4], it[5] ], [ it[0], it[6], it[7] ] ] } //rearange for convenience
      .combine( bam2match_ip_mut.map{ it-> [ it[0], it[1],[ it[2][1], it[2][2] ] ] }, by:0 )  //match bam with IP_mut
      .combine( bam2match_input_mut.map{ it-> [ it[1], it[0],[ it[2][1], it[2][2] ] ] }, by:1 ) //match bam with INPUT_mit
      .map{ it  -> [ it[2], it[3], [ it[1], it[4], it[5] ], [it[0], it[6], it[7] ] ] } //rearange for convenience -> [ [IP_w], [INPUT_w], [IP_m], [INPUT_m] ]
      .into{ bam_WT_vs_MUT_2bamIP ; bam_WT_vs_MUT }



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       BAM INDEXING                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


bam_WT_vs_MUT_2bamIP
        .flatMap{ sample -> [ sample[0], sample[2] ] }
        .unique()
        .into{ bam_IP_2index; bam_IP_2joinBAI }



process IP_bam_indexing {
  tag "$file_id"
  label "sambamba"
  input:
    set file_id, data_type , file(bams) from bam_IP_2index

  output:
    set file_id, "*.bam*" into bai_2join_bam_IP

  script:
  """
  sambamba index -t ${task.cpus} ${bams[0]}
  sambamba index -t ${task.cpus} ${bams[1]}
  """
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   CONVERT BAM 2 BIGWIG                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


bam_IP_2joinBAI.join(bai_2join_bam_IP, by: 0).set{ bam_bai_IP_2bigwig }


process IP_bam2BigWig {
  tag "$file_id"
  label "deeptools"

  publishDir "${params.outdir}/bigWig/${file_id}/", mode: 'copy'

  input:
    set file_id , data_type , file(bams), file(idx) from bam_bai_IP_2bigwig


  output:
    set file_id, "*_fasta_occupancy.bw", '*_fastaCalib_occupancy.bw' into ip_bw_2norm, ip_bw_2stats

  script:
    """
    bamCoverage --binSize 1 -p ${task.cpus} --ignoreDuplicates -b ${bams[0]} -o ${file_id}_fasta_occupancy.bw
    bamCoverage --binSize 1 -p ${task.cpus} --ignoreDuplicates -b ${bams[1]} -o ${file_id}_fastaCalib_occupancy.bw
    """
}



if (! params.duplicate_removal){
    mapping_reports2normalization.into{ log_2match_ip_wt ; log_2match_input_wt; log_2match_ip_mut ; log_2match_input_mut}

    experiment_ID_2log_2IP_vs_INPUT
          .combine( log_2match_ip_wt.map{ it-> [ it[0], it[1], it[3] ] }, by:0 )  //match log with IP_wt
          .combine( log_2match_input_wt.map{ it-> [ it[1], it[0], it[3] ] }, by:1 ) //match log with INPUT_wt
          .map{ it  -> [ it[2], it[3],  it[1], [ it[4], it[5] ] ,  it[0], [ it[6], it[7] ]  ] } //rearange for convenience
          .combine( log_2match_ip_mut.map{ it-> [ it[0], it[1], it[3] ] }, by:0 )  //match log with IP_mut
          .combine( log_2match_input_mut.map{ it-> [ it[1], it[0], it[3] ] }, by:1 )  //match log with INPUT_mit
          .map{ it  -> [ it[2], it[3],  it[4],it[5],  it[1], [ it[6], it[7] ] , it[0], [ it[8], it[9] ]  ] } //rearange for convenience -> [ [IP_w], [INPUT_w], [IP_m], [INPUT_m] ]
          .flatMap{ sample -> [ [ sample[0], sample[1], sample[2], sample[3] ], [ sample[4], sample[5], sample[6], sample[7] ] ] }
          .unique()
          .set{ log_2IP_vs_INPUT }


}
else{

  dedup_reports2normalization.into{ log_2match_ip_wt ; log_2match_input_wt; log_2match_ip_mut ; log_2match_input_mut}

  experiment_ID_2log_2IP_vs_INPUT
        .combine( log_2match_ip_wt.map{ it-> [ it[0], it[1], it[3] ] }, by:0 )  //match log with IP_wt
        .combine( log_2match_input_wt.map{ it-> [ it[1], it[0], it[3] ] }, by:1 ) //match log with INPUT_wt
        .map{ it  -> [ it[2], it[3],  it[1], [ it[4], it[5] ] ,  it[0], [ it[6], it[7] ]  ] } //rearange for convenience
        .combine( log_2match_ip_mut.map{ it-> [ it[0], it[1], it[3] ] }, by:0 )  //match log with IP_mut
        .combine( log_2match_input_mut.map{ it-> [ it[1], it[0], it[3] ] }, by:1 )  //match log with INPUT_mit
        .map{ it  -> [ it[2], it[3],  it[4],it[5],  it[1], [ it[6], it[7] ] , it[0], [ it[8], it[9] ]  ] } //rearange for convenience -> [ [IP_w], [INPUT_w], [IP_m], [INPUT_m] ]
        .flatMap{ sample -> [ [ sample[0], sample[1], sample[2], sample[3] ], [ sample[4], sample[5], sample[6], sample[7] ] ] }
        .unique()
        .set{ log_2IP_vs_INPUT}


}


//v.view()

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   NORMALIZATION                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

ip_bw_2norm
          .join(log_2IP_vs_INPUT)
          .set{ input_norm_bw }


process Normalized_bigWig {
  label 'R'
 tag "$file_id_IP"
 publishDir "${params.outdir}/bigWig/${file_id_IP}/norm", mode: 'copy'

 input:
    set file_id_IP, file(ip_fasta_occupancy_bw), file(ip_fastaCalib_occupancy_bw), file(log_IP), file_id_INPUT, file(log_INPUT) from input_norm_bw


 output:
    set file_id_IP, "*_fasta_occupancy_norm.bw", "*_fastaCalib_occupancy_norm.bw" into norm_ip_bw_2stats

 script:
 if (params.shortReads){
   tool=1
 }
 if (params.duplicate_removal){
   tool="RD"
}
if(!params.duplicate_removal && !params.shortReads){
   tool=2
 }
   """
   get_normalization_factor.R \
   ${log_IP[0]} \
   ${log_IP[1]} \
   ${log_INPUT[0]} \
   ${log_INPUT[1]} ${tool} > norm_factor

   normalized_occupancy.R \
   ${ip_fasta_occupancy_bw} \
   ${file_id_IP}_fasta_occupancy_norm.bw \
   norm_factor

   normalized_occupancy.R \
   ${ip_fastaCalib_occupancy_bw} \
   ${file_id_IP}_fastaCalib_occupancy_norm.bw \
   norm_factor
   """
 }


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --               NORMALIZATION STATISTICS                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////




if (params.skipNorm_stats){
  Channel
        .empty()
        .set { plots_occupancy }

}

else{

  ip_bw_2stats
        .into{ ip_bw_2match_IP_w; ip_bw_2match_IP_m }

  norm_ip_bw_2stats
        .into{ norm_IP_bw_2match_IP_w; norm_IP_bw_2match_IP_m }

  experiment_ID_2bw_IP_WT_MUT
              .combine( ip_bw_2match_IP_w.map{ it-> [ it[0], it[1], it[2] ] }, by:0 )  //match bw  with IP_wt
              .combine( ip_bw_2match_IP_m.map{ it-> [ it[1], it[2], it[0] ] }, by:2 ) //match  bw with IP_mut
              .map{ it  -> [ it[1], it[5], it[0], it[7] ] } //rearange for convenience
              .set{ bw_IP_C_WT_MUT }


  experiment_ID_2bw_norm_IP_WT_MUT
              .combine( norm_IP_bw_2match_IP_w.map{ it-> [ it[0], it[1], it[2] ] }, by:0 )  //match bw  with IP_wt
              .combine( norm_IP_bw_2match_IP_m.map{ it-> [ it[1], it[2], it[0] ] }, by:2 ) //match  bw with IP_mut
              .map{ it  -> [ it[1], it[5], it[0], it[7] ] } //rearange for convenience
              .set{ bw_norm_IP_C_WT_MUT }

  bw_IP_C_WT_MUT.concat(bw_norm_IP_C_WT_MUT).set{ ip_C_WT_MUT_bw }

process Statistics_normalization {
  label 'R'
 tag "$file_id_ip_c_wt"
 publishDir "${params.outdir}/bigWig/statistics/${file_id_ip_c_wt}_vs_${file_id_ip_c_mut}", mode: 'copy'


 input:
   set file_id_ip_c_wt, file(wt_occupancy_bw), file_id_ip_c_mut, file(mut_occupancy_bw) from ip_C_WT_MUT_bw

 output:
    file "*.svg" into plots_occupancy

 script:
   """
   occupancy_statistics.R \
   ${wt_occupancy_bw} \
   ${mut_occupancy_bw} \
   ${wt_occupancy_bw.baseName}_vs_${mut_occupancy_bw.baseName}.svg
   """
}
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     MERGING BAM                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



process Merging_BAM {

 tag "$ip_id_merged"
 publishDir "${params.outdir}/mapping/BAM_merged/${ip_id_merged}", mode: 'copy'
 label "samtools"

 input:
   set ip_wt, input_wt, ip_mut, input_mut from bam_WT_vs_MUT

 output:
  set ip_id_merged , input_id_merge, file("${file_id_IP_wt}_${file_id_IP_mut}.bam"), file("${file_id_INPUT_wt}_${file_id_INPUT_mut}.bam"), val(data_type) into bam_IP_vs_INPUT_merged

 script:
 file_id_IP_wt=ip_wt[0]
 data_type=ip_wt[1]
 bam_ip_wt=ip_wt[2][0]

 file_id_INPUT_wt=input_wt[0]
 bam_input_wt=input_wt[2][0]

 file_id_IP_mut=ip_mut[0]
 bam_ip_mut=ip_mut[2][0]

 file_id_INPUT_mut=input_mut[0]
 bam_input_mut=input_mut[2][0]

 ip_id_merged = "${file_id_IP_wt}_${file_id_IP_mut}"
 input_id_merge = "${file_id_INPUT_wt}_${file_id_INPUT_mut}"

   """
   samtools merge ${file_id_IP_wt}_${file_id_IP_mut}.bam ${bam_ip_wt} ${bam_ip_mut}
   samtools merge ${file_id_INPUT_wt}_${file_id_INPUT_mut}.bam ${bam_input_wt} ${bam_input_mut}
   """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       PEAK CALLING                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Peak_calling {
  tag "${ip_ID} vs ${input_ID}"
  label "macs2"
  publishDir "${params.outdir}/peakCalling/${ip_ID}", mode: 'copy'

  input:
    set ip_ID, input_ID, file(bam_ip), file(bam_input), data_type from bam_IP_vs_INPUT_merged

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
     file report_ddup from dDup_bam_report.collect().ifEmpty([])
     file report_peakC from peak_calling_report.collect().ifEmpty([])

     output:
     file "*multiqc_*" into multiqc_report


     when:
     !params.skipMultiqc

     script:
     """
     multiqc -f . \\
     -m fastqc -m cutadapt -m picard -m macs2
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
    ${c_purple}  QUANTITATIVE CHIP-SEQ Pipeline          ${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
