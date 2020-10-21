log.info "fastq files : ${params.fastq}"
log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_files }
Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .set { bed_files }
Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process adaptor_removal {
  tag "$pair_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, "*_cut_R{1,2}.fastq.gz" into fastq_files_cut

  script:
  """

  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
  -o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${pair_id}_report.txt
  """
}

process trimming {
  tag "${reads}"
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_cut

  output:
  set pair_id, "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
> ${pair_id}_trimming_report.txt
"""
}

process fasta_from_bed {
  tag "${bed.baseName}"
  publishDir "results/fasta/", mode: 'copy'

  input:
  file fasta from fasta_files
  file bed from bed_files

  output:
  file "*_extracted.fasta" into fasta_files_extracted

  script:
"""
bedtools getfasta -name \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}

process index_fasta {
  tag "$fasta.baseName"
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_files_extracted

  output:
    file "*.index*" into index_files
    file "*_kallisto_report.txt" into index_files_report

  script:
"""
kallisto index -k 31 --make-unique -i ${fasta.baseName}.index ${fasta} \
2> ${fasta.baseName}_kallisto_report.txt
"""
}

process mapping_fastq {
  tag "$reads"
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_trim
  file index from index_files.collect()

  output:
  file "*" into counts_files

  script:
"""
mkdir ${pair_id}

kallisto quant -i ${index} -t ${task.cpus} \
--bias --bootstrap-samples 100 -o ${pair_id} \
${reads[0]} ${reads[1]} &> ${pair_id}/kallisto_report.txt
"""
}
