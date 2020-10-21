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
      --inputFromSRA                Specify SRA IDs fastq files data (required if --fromPATH not specified)
      --inputFromPATH               Directory pattern for fastq files: (required if --fromSRA not specified)
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
      --quantile_norm               Normalized each position by quantile (not yet got it right) (default: average occupancy normalization)

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
/* --                SETUP CONFIGURATION VARIABLES                       -- */
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
params.inputFromPATH = false
params.quantile_norm = false
params.inputFromSRA = false
params.fasta = false
params.outdir = 'results'


/*
 * SETUP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

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
else { exit 1,
  log.warn "=================================================================\n" +
           "  WARNING! No genome fasta file precised.\n" +
           "  Use '--fasta' \n" +
           "  Or '--help' for more informations"
           "======================================================================="
}


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
            .fromFilePairs(params.inputFromPATH, size:-1)
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
  summary['Normalization']          = params.quantile_norm ? 'quantile' : 'average occupancy'
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
        file "${fasta.simpleName}.fasta" into fasta_2indexing, fasta_2sizeChr

      script:
      """
      zcat ${fasta} > ${fasta.simpleName}.fasta
      """
     }
 }
 else {
   fasta_file_zip.into { fasta_2sizeChr;
                       fasta_2indexing }
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
    tag "$file_id"

    if (!params.shortReads){
      label "bowtie2"
    }
    else {
      label "bowtie"
    }

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
/* --     GET CHR SIZE FOR DANPOS QUANTILE NORMALIZATION                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



process Get_size_chromosome {
     tag "$fasta.baseName"
     publishDir "${params.outdir}/genome/", mode: 'copy'
     label "samtools"
     input:
       file fasta from fasta_2sizeChr

     output:
       file "sizes.genome" into chr_size_file2norm, chr_size_file2bw

     script:
       """
           samtools faidx ${fasta}
           cut -f1,2 ${fasta}.fai > sizes.genome
       """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         PEAK CALLING                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Peak_calling {
  label "danpos"
  publishDir "${params.outdir}/peakCalling/", mode: 'copy'
  tag "$file_id"

  input:
  set file_id, data_type, file(bam) from dedup_bam_files


  output:
  set file_id, "${file_id}/*.wig" into wig_files, wig_files_2concat
  set file_id, "${file_id}/*.xls" into xls_files

  script:
  m = data_type=="SE" ? 0 : 1
  """
  danpos.py dpos -o "./${file_id}" -a 1 -n N -z 0 --paired $m ${bam}
  mv ./${file_id}/pooled/*.xls ./${file_id}/
  mv ./${file_id}/pooled/*.wig ./${file_id}/
  rm -r ./${file_id}/pooled/
  """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     NORMALIZED WIG FILE                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Wig_normalization {
  if (! params.quantile_norm){
    label "wig_manipulation"
  }
  else {
    label "danpos"
  }

  publishDir "${params.outdir}/peakCalling/${file_id}/normalization", mode: 'copy'

  tag "$file_id"

  input:
  set file_id, file(wig) from wig_files
  file(chr_size) from chr_size_file2norm.collect().ifEmpty([])

  output:
  set file_id, "*nor.wig" into wig_norm_files, wig_norm_files_2concat

  script:
  if (! params.quantile_norm){
  """
  wig_average_normalization.py ${wig} ${file_id}.anor.wig
  """
  }
  else{
  """
  danpos.py wiq --out_dir  "./norm" --buffer_size 50 --step 1 ${chr_size} ${wig}
  mv ./norm/* .
  """
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         SMOOTH WIG FILE                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Wig_smooth {
  label "danpos"
  publishDir "${params.outdir}/peakCalling/${file_id}/normalization", mode: 'copy'

  tag "$file_id"

  input:
  set file_id, file(wig) from wig_norm_files


  output:
  set file_id, "*nor.smooth.wig" into wig_norm_smooth_files, wig_norm_smooth_files_2concat
  set file_id, "*nor.smooth.positions.xls" into xls_norm_files

  script:
  """
  danpos.py dpos -a 1 -n N ${wig}
  mv result/pooled/${file_id}*.wig .
  mv result/pooled/${file_id}*.xls .
  rm -r ./result/
    """
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --               BUILD BED FILE WITH NDR POSITIONS                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
 * Here we defined nucleosome-depleted regions
 * (NDRs) as regions spanning at least 150 nucleotides (corresponding
 * to the eviction of at least one nucleosome) with a normalized sequence
 * coverage lower than 0.4.
 */

process NDRs_finder {
  label "wig_manipulation"
  publishDir "${params.outdir}/NDRs/${file_id}", mode: 'copy'
  tag "$file_id"

  input:
  set file_id, file(wig_norm) from wig_norm_smooth_files

  output:

  set file_id, "*.bed" into bed_files

  script:
  if (! params.quantile_norm){
  """
  wig2bed.py ${wig_norm} 150 0.4 ${file_id}.bed
  """
  }
  else{
  """
  wig2bed.py ${wig_norm} 150 0.4 ${file_id}.bed
  """
  }
}

wig_files_2concat
        .concat(wig_norm_files_2concat)
        .concat(wig_norm_smooth_files_2concat)
        .set { wig_files_2bw }


process Wig2bigWig {
  label "wig2bigwig"
  publishDir "${params.outdir}/BigWig/${file_id}", mode: 'copy'
  tag "$file_id"

  input:
  set file_id, file(wig) from wig_files_2bw
  file(chr_size) from chr_size_file2bw.collect()


  output:
  set file_id, "*.bw" into bw_files

  script:
  """
  wigToBigWig -clip ${wig} ${chr_size} ${wig.baseName}.bw
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

     output:
     file "*multiqc_*" into multiqc_report


     when:
     !params.skipMultiqc

     script:
     """
     multiqc -f . \\
     -m fastqc -m cutadapt -m bowtie1 -m bowtie2 -m picard
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
    ${c_purple}  MNASE-SEQ Pipeline          ${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
