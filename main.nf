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

// Validate inputs
Channel.fromPath(params.batch)
    .ifEmpty { exit 1, "Batch CSV file not found: ${params.batch}" }
    .splitCsv(skip: 1)
    .map { molecular_id, omics_sample_name, dir, key -> [file(key).baseName, file(key), molecular_id, omics_sample_name, dir] }
    .set { batch }
batch
    .map { name, bam, molecular_id, omics_sample_name, dir -> [name, bam] }
    .into { batch_sam_to_fastq; batch_merge_bams }
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_merge_bams }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
           .set { fai_merge_bams }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .set { dict_merge_bams }
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
  set val(name), file(bam) from batch_sam_to_fastq

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

bwa_index = fasta_bwa.merge(bwa_index_amb,bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = reads.combine(bwa_index)

/*--------------------------------------------------
  Align reads to reference genome
---------------------------------------------------*/

process BwaMem {
  tag "${name}"

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
  memory "14.GB"
  cpus 4

  input:
  set val(name), file(reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

  output:
  set val(name), file("${name}.aligned.bam") into aligned_bam

  script:
  """
  /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${fasta} ${reads} | samtools view -1bS > ${name}.aligned.bam
  """
}

aligned_bam_ref = aligned_bam.merge(fasta_merge_bams, fai_merge_bams, dict_merge_bams)
batch_merge_bams.combine(aligned_bam_ref, by: 0)
  .set { merge_bams }

/*--------------------------------------------------
  Merge original input uBAM with BWA-aligned BAM
---------------------------------------------------*/

process MergeBamAlignment {
  tag "${name}"
  publishDir "${params.outdir}", mode: 'copy'

  container 'broadinstitute/gatk:4.0.4.0'
  memory "8.GB"

  input:
  set val(name), file(unmapped_bam), file(aligned_bam), file(fasta), file(fai), file(dict) from merge_bams

  output:
  file("${name}.merged.bam") into merged_bam

  script:
  """
  /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      MergeBamAlignment \
     --ALIGNED_BAM ${aligned_bam} \
     --UNMAPPED_BAM ${unmapped_bam} \
     --OUTPUT ${name}.merged.bam \
     --REFERENCE_SEQUENCE ${fasta} \
     --PAIRED_RUN true \
     --SORT_ORDER coordinate \
     --CREATE_INDEX true \
     --CLIP_ADAPTERS true \
     --MAX_RECORDS_IN_RAM 2000000 \
     --MAX_INSERTIONS_OR_DELETIONS -1 \
     --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
     --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
     --ALIGNER_PROPER_PAIR_FLAGS true
  """
}