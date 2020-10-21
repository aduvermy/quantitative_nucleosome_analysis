#!/usr/bin/env nextflow
/*
========================================================================================
                        mnase_analysis
========================================================================================
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run src/mnase_analysis.nf -c src/mnase_analysis.config --fromSRA ["A","B","C"] --fromPATH ["A","B","C"] --fasta <genome file> -profile singularity


    Required arguments:
      --inputFromSRA                     Specify SRA IDs fastq files data (required if --fromPATH not specified)
      --inputFromPATH                    Directory pattern for fastq files: (required if --fromSRA not specified)
                                    Please use following syntax for paired end data : data/filePE_R{1,2}*.fastq
    Reference genome
      --fasta                       Full path to directory containing genome fasta file

    Mapping option:
      --shortReads                  Specifies that all input file are not long reads in order to use more adapted mapper: bowtie (default: bowtie2)

    QC Option:
      --skipFastqc                  Skip reads quality control step (default: activated).
      --skipMultiqc                 Skip merging tools reports suitable with multiqc (default: activated)

    Trimming option:
      --trimming                    Activated trimming step (default: desactivated)
      --adapter_removal             Activated adapter removal step (default: desactivated)

    Remove duplicates
      --duplicate_removal           Activated reads duplicates removal step (default: desactivated)

    Peak Calling normalization:
      --average_occ_norm            Normalized each position by average occupancy (default: quantile normalization)

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

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

//params.inputFromSRA = ["SRR7289789", "SRR7289788", "SRR7289787", "SRR7289786"]
params.inputFromPATH=['data/simulation/l50/mix25_75.fastq', 'data/simulation/l50/mix25_75_mut.fastq', 'data/simulation/l50/mix40_60_mut.fastq', 'data/simulation/l50/mix50_50.fastq']
params.fasta_calib="/Xnfs/lbmcdb/common/Genomes/Schizosaccharomyces_cerevisiae/S288C_reference_sequence_R64-2-1_20150113.fsa"
params.fasta="/Xnfs/lbmcdb/common/Genomes/Schizosaccharomyces_pombe/2017_09_19_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz"
params.outdir = './results-mapping_investigation_simul'



////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

if (!params.inputFromSRA && !params.inputFromPATH){  exit 1, "Input params not specified!\nUse --help" }

if (params.fasta) {
        Channel
            .fromPath( params.fasta )
            .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
            .set { fasta_file }
}
else { exit 1, "FASTA genome file not specified!" }

if (params.fasta_calib) {
        Channel
            .fromPath( params.fasta_calib )
            .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
            .set { fasta_calib_file }
}
else { exit 1, "FASTA genome file for calibration not specified!" }




if (params.inputFromSRA){
        Channel
            .fromSRA(params.inputFromSRA, apiKey:'6e15df3377f722be16ef0e546d8a40982808')
            .ifEmpty { error "Cannot find any SRA IDs matching: ${params.inputFromSRA}" }
            .set{ fastq_fromSRA }
}
else{
        Channel
            .empty()
            .set { fastq_fromSRA }

}

if (params.inputFromPATH){
        Channel
            .fromFilePairs(params.inputFromPATH, size: -1)
            .ifEmpty { error "Cannot find any file matching: ${params.inputFromPATH}" }
            .set{ fastq_fromPATH }
}
else{
        Channel
            .empty()
            .set { fastq_fromPATH }
}


/*
* CONCAT fromPath and fromSRA channels into a single input channel
*/
fastq_fromPATH
              .concat(fastq_fromSRA)
              .into{ fastq_raw_2QC ; fastq_raw_2trim }




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
  summary['From SRA']               = params.inputFromSRA ? params.inputFromSRA : 'Not supplied'
  summary['From PATH']              = params.inputFromPATH ? params.inputFromPATH : 'Not supplied'
  summary['Remove Adapter']         = params.adapter_removal ? 'Yes' : 'Skipped'
  summary['Trimming']               = params.trimming ? 'Yes': 'Skipped'
  summary['Reads QC']               = params.skipFastqc ? 'Skipped' : 'Yes'
  summary['Merging Reports']        = params.skipMultiqc ? 'Skipped' : 'Yes'
  summary['Mapper']                 = params.shortReads ? 'Bowtie1' : 'Bowtie2'
  summary['Remove Duplicate']       = params.duplicate_removal ? 'Yes': 'Skipped'
  summary['Config Profile']         = workflow.profile
  summary['Output']                 = params.outdir
  log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
  log.info "-\033[2m--------------------------------------------------\033[0m-"


  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  /* --                                                                     -- */
  /* --                         TRIMMING READS                              -- */
  /* --                                                                     -- */
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  if (!params.trimming) {
       fastq_raw_2trim.set{fastq_trim_2cutAdapt }
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
   fasta_file_zip.set{ fasta_2indexing }
 }


 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 /* --                                                                     -- */
 /* --                   UNZIP FASTA CALIB FILE IF NEEDED                        -- */
 /* --                                                                     -- */
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////

 /*
  * GET EXTENSION FASTA FILE TO TEST IT
  */

 fasta_calib_file.into{ fasta_calib_file_test_zip;
                  fasta_calib_file_zip }

 ext=fasta_calib_file_test_zip.getVal().getExtension()


  if(ext=="gz" || ext=="bz" || ext=="zip"){

     log.info "Genome fasta calib file zip, we need to unzip it for next steps"

     process Unzip_fasta_calib {
      tag "$fasta.simpleName"

      input:
        file fasta from fasta_calib_file_zip

      output:
        file "${fasta.simpleName}.fasta" into fasta_calib_2indexing

      script:
      """
      zcat ${fasta} > ${fasta.simpleName}.fasta
      """
     }
 }
 else {
   fasta_calib_file_zip.set{ fasta_calib_2indexing }
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
 /* --                      BUILD GENOME INDEX                             -- */
 /* --                                                                     -- */
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////

 process Fasta_calib_indexing {
     if (!params.shortReads){
           label "bowtie2"
     }
     else {
           label "bowtie"
     }
     tag "$fasta.simpleName"

     input:
     file fasta from fasta_calib_2indexing

     output:
     file "*.index*" into index_calib_files
     file "*_report.txt" into indexing_calib_report

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
    file index from index_files.collect()


    output:
    set file_id, val(data_type) ,"*.bam" into bam_files
    set file_id, "*.fastq" into reads_unaligned_2fasta
    //set file_id, "*_aligned_2fasta.fastq" into reads_aligned_2fasta
    file "*.txt" into mapping_report

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
        """
      }
      else{
        data_type="PE"
        """
        # -v specify the max number of missmatch, -k the number of match reported per
        # reads
        bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
        --un ${file_id}_unaligned_2fasta_calib.fastq --al ${file_id}_aligned_2fasta_calib.fastq \
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


        """
      }
      else{
        data_type="PE"
        """
        bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
        --un ${file_id}_unaligned_2fasta.fastq --al ${file_id}_aligned_2fasta.fastq \
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

//bedtools genomecov -ibam SRR7289786_exclusiv2fasta.bam -g /scratch/Bio/aduvermy/quantitative-nucleosome-analysis/data/genome/S.pombesizes.genome -d
//picard-tools MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="Saduvermy@m6142comp2:/scratch/Bio/aduvermy/quantitative-nucleosome-analysis$ picard-tools MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="SRR7289789_exclusiv2fastaCalib.bam" OUTPUT="SRR7289789_exclusiv2fastaCalib_dedup.bam" METRICS_FILE=SRR7289789_exclusiv2fastaCalib_metrics.txt REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE PROGRAM_RECORD_ID=null


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
    ${c_purple}  MAPPING Pipeline          ${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
