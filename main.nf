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

/*--------------------------------------------------
  Read unmapped BAM, convert to FASTQ
---------------------------------------------------*/

process SamToFastq {
  tag "${name}"
  publishDir "${params.outdir}", mode: 'copy'

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'

  input:
  set val(molecular_id), val(omics_sample_name), val(dir), val(name), file(bam) from batch

  output:
  file("${name}.fastq") into reads

  script:
  """
  java -Dsamjdk.compression_level=5 -Xms4g -jar /usr/gitc/picard.jar \
      SamToFastq \
      I=${bam} \
      FASTQ=${name}.fastq \
      INTERLEAVE=true
  """
}