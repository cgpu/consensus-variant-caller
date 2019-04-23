#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/consensus-variant-caller
========================================================================================
 lifebit-ai/consensus-variant-caller Consensus variant calling workflow for human panel-based DNA sequencing
 Input requirements:
 - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
   - filenames all have the same suffix (we use ".unmapped.bam")
   - files must pass validation by ValidateSamFile
   - reads are provided in query-sorted order
   - all reads must have an RG tag

 Output :
 - recalibrated bam and it's index and md5
 - GATK vcf
 - samtools/bcftools vcf
 - Annovar annotated vcfs and tabular file for each variant caller

 Software version requirements (see recommended dockers in inputs JSON)
 - GATK 4 or later (see gatk docker)
 - Picard (see gotc docker)
 - Samtools (see gotc docker)

 #### Homepage / Documentation
 https://github.com/lifebit-ai/fred-hutch-gatk
----------------------------------------------------------------------------------------
*/

Channel.fromPath(params.batch)
    .ifEmpty { exit 1, "Batch CSV file not found: ${params.batch}" }
    .splitCsv(skip: 1)
    .map { molecular_id, omics_sample_name, dir, key -> [molecular_id, omics_sample_name, dir, file(key).baseName, file(key)] }
    .set { batch }
// Validate inputs
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .set { fasta_bwa }
}
params.bwa_index_amb = params.genome ? params.genomes[ params.genome ].bwa_index_amb ?: false : false
if (params.bwa_index_amb) {
    Channel.fromPath(params.bwa_index_amb)
           .ifEmpty { exit 1, "bwa_index_amb annotation file not found: ${params.bwa_index_amb}" }
           .set { bwa_index_amb }
}
params.bwa_index_ann = params.genome ? params.genomes[ params.genome ].bwa_index_ann ?: false : false
if (params.bwa_index_ann) {
    Channel.fromPath(params.bwa_index_ann)
           .ifEmpty { exit 1, "bwa_index_ann annotation file not found: ${params.bwa_index_ann}" }
           .set { bwa_index_ann }
}
params.bwa_index_bwt = params.genome ? params.genomes[ params.genome ].bwa_index_bwt ?: false : false
if (params.bwa_index_bwt) {
    Channel.fromPath(params.bwa_index_bwt)
           .ifEmpty { exit 1, "bwa_index_bwt annotation file not found: ${params.bwa_index_bwt}" }
           .set { bwa_index_bwt }
}
params.bwa_index_pac = params.genome ? params.genomes[ params.genome ].bwa_index_pac ?: false : false
if (params.bwa_index_pac) {
    Channel.fromPath(params.bwa_index_pac)
           .ifEmpty { exit 1, "bwa_index_pac annotation file not found: ${params.bwa_index_pac}" }
           .set { bwa_index_pac }
}
params.bwa_index_sa = params.genome ? params.genomes[ params.genome ].bwa_index_sa ?: false : false
if (params.bwa_index_sa) {
    Channel.fromPath(params.bwa_index_sa)
           .ifEmpty { exit 1, "bwa_index_sa annotation file not found: ${params.bwa_index_sa}" }
           .set { bwa_index_sa }
}

/*--------------------------------------------------
  Read unmapped BAM, convert to FASTQ
---------------------------------------------------*/

process SamToFastq {
  tag "${name}"

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
  memory "14.GB"
  cpus 4

  input:
  set val(molecular_id), val(omics_sample_name), val(dir), val(name), file(bam) from batch

  output:
  set val(name), file("${name}.fastq") into reads

  script:
  """
  java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
      INPUT=${bam} \
      FASTQ=${name}.fastq \
      INTERLEAVE=true \
      NON_PF=true 
  """
}

bwa = reads.merge(fasta_bwa, bwa_index_amb, bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)

/*--------------------------------------------------
  Align reads to reference genome
---------------------------------------------------*/

process BwaMem {
  tag "${name}"
  publishDir "${params.outdir}", mode: 'copy'

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
  memory "14.GB"
  cpus 4

  input:
  set val(name), file(reads), file(fasta), file(bwa_index) from bwa

  output:
  file("${name}.aligned.bam") into aligned_bam

  script:
  """
  /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${fasta} ${reads} | samtools view -1bS > ${name}.aligned.bam
  """
}