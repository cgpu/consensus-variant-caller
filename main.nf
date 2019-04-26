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
           .into { fasta_bwa; fasta_ref }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
           .set { fai_ref }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .set { dict_ref }
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
params.dbsnp_gz = params.genome ? params.genomes[ params.genome ].dbsnp_gz ?: false : false
if (params.dbsnp_gz) {
    Channel.fromPath(params.dbsnp_gz)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
           .set { dbsnp_gz}
}
params.dbsnp_idx_gz = params.genome ? params.genomes[ params.genome ].dbsnp_idx_gz ?: false : false
if (params.dbsnp_idx_gz) {
    Channel.fromPath(params.dbsnp_idx_gz)
           .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
           .set { dbsnp_idx_gz}
}
params.golden_indel_gz = params.genome ? params.genomes[ params.genome ].golden_indel_gz ?: false : false
if (params.golden_indel_gz) {
    Channel.fromPath(params.golden_indel_gz)
           .ifEmpty { exit 1, "golden_indel_gz annotation file not found: ${params.golden_indel_gz}" }
           .set { golden_indel_gz }
}
params.golden_indel_idx_gz = params.genome ? params.genomes[ params.genome ].golden_indel_idx_gz ?: false : false
if (params.golden_indel_idx_gz) {
    Channel.fromPath(params.golden_indel_idx_gz)
           .ifEmpty { exit 1, "golden_indel_idx_gz annotation file not found: ${params.golden_indel_idx_gz}" }
           .set { golden_indel_idx_gz }
}
Channel.fromPath(params.annotation)
    .ifEmpty { exit 1, "Annotation folder not found: ${params.annotation}" }
    .set { annotation }

fasta_ref.merge(fai_ref, dict_ref)
  .into { merge_bams_ref; baserecalibrator_ref; applybqsr_ref; haplotypecaller_ref; bcftools_ref}

/*--------------------------------------------------
  Gunzip the dbSNP VCF file if gzipped
---------------------------------------------------*/

process gunzipDbsnp {
    tag "$dbsnp_gz"

    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'

    input:
    file dbsnp_gz from dbsnp_gz
    file dbsnp_idx_gz from dbsnp_idx_gz

    output:
    file "*.vcf" into dbsnp
    file "*.vcf.idx" into dbsnp_idx

    script:
    if ( "${dbsnp_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $dbsnp_gz
     gunzip -d --force $dbsnp_idx_gz
     """
   } else {
     """
     cp $dbsnp_gz dbsnp.vcf
     cp $dbsnp_idx_gz dbsnp.vcf.idx
     """
   }
}

/*--------------------------------------------------
  Gunzip the golden indel file if gzipped
---------------------------------------------------*/

process gunzipGoldenIndel {
  tag "$golden_indel_gz"

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'

  input:
  file golden_indel_gz from golden_indel_gz
  file golden_indel_idx_gz from golden_indel_idx_gz

  output:
  file "*.vcf" into golden_indel
  file "*.vcf.idx" into golden_indel_idx

  script:
  if ( "${golden_indel_gz}".endsWith(".gz") ) {
    """
    gunzip -d --force $golden_indel_gz
    gunzip -d --force $golden_indel_idx_gz
    """
  } else {
    """
    cp $golden_indel_gz golden_indel.vcf
    cp $golden_indel_idx_gz golden_indel.vcf.idx
    """
  }
}

/*--------------------------------------------------
  Read unmapped BAM, convert to FASTQ
---------------------------------------------------*/

process samToFastq {
  tag "${name}"

  container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
  memory "14.GB"
  cpus 4

  input:
  set val(name), file(bam) from batch_sam_to_fastq

  output:
  set val(name), file("${name}.fastq") into reads, reads_fastqc

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

/*--------------------------------------------------
  Use FASTQC for FASTQ quality control
---------------------------------------------------*/

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
    container 'flowcraft/fastqc:0.11.7-1'

    input:
    set val(name), file(reads) from reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    when: !params.skip_multiqc

    script:
    """
    fastqc -q $reads
    """
}

bwa_index = fasta_bwa.merge(bwa_index_amb,bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = reads.combine(bwa_index)

/*--------------------------------------------------
  Align reads to reference genome
---------------------------------------------------*/

process bwaMem {
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

bams_to_merge = batch_merge_bams.combine(aligned_bam, by: 0)
merge_bams = bams_to_merge.combine(merge_bams_ref)

/*--------------------------------------------------
  Merge original input uBAM with BWA-aligned BAM
---------------------------------------------------*/

process mergeBamAlignment {
  tag "${name}"

  container 'broadinstitute/gatk:4.0.4.0'
  memory "8.GB"

  input:
  set val(name), file(unmapped_bam), file(aligned_bam), file(fasta), file(fai), file(dict) from merge_bams

  output:
  set val(name), file("${name}.merged.bam"), file("${name}.merged.bai") into merged_bam, merged_bam_applybqsr, merged_bam_qc

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

/*--------------------------------------------------
  Run BAM QC using Qualimap before recalibration
---------------------------------------------------*/

process runBamQCmapped {
    tag "$bam"
    container 'maxulysse/sarek:latest'

    input:
    set val(name), file(bam), file(index) from merged_bam_qc

    output:
    file("${name}_mapped") into bamqc_mapped_report

    when: !params.skip_multiqc

    script:
    """
    qualimap \
      bamqc \
      -bam ${bam} \
      --paint-chromosome-limits \
      --genome-gc-distr HUMAN \
      -nt ${task.cpus} \
      -skip-duplicated \
      --skip-dup-mode 0 \
      -outdir ${name}_mapped \
      -outformat HTML
    """
}

baserecalibrator_ref = baserecalibrator_ref.merge(dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = merged_bam.combine(baserecalibrator_ref)

/*--------------------------------------------------
  Generate Base Quality Score Recalibration (BQSR) model
---------------------------------------------------*/

process baseRecalibrator {
  tag "${name}"

  container 'broadinstitute/gatk:4.0.4.0'
  memory "8G"

  input:
  set val(name), file(bam), file(bai), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

  output:
  set val(name), file("${name}.recal_data.table") into baserecalibrator_table
  
  script:
  // TODO: add BED file?
  """
  /gatk/gatk --java-options "-Xms4g" \
      BaseRecalibrator \
      -R ${fasta} \
      -I ${bam} \
      --use-original-qualities \
      -O ${name}.recal_data.table \
      --known-sites ${dbsnp} \
      --known-sites ${golden_indel}
  """
}

applybqsr_bam_table = merged_bam_applybqsr.combine(baserecalibrator_table, by: 0)
applybqsr = applybqsr_bam_table.combine(applybqsr_ref)

/*--------------------------------------------------
  Apply Base Quality Score Recalibration (BQSR) model
---------------------------------------------------*/

process applyBQSR {
  tag "${name}"

  container 'broadinstitute/gatk:4.0.4.0'
  memory "8G"

  input:
  set val(name), file(bam), file(bai), file(recalibration_table), file(fasta), file(fai), file(dict) from applybqsr

  output:
  set val(name), file("${name}.recal.bam"), file("${name}.recal.bai") into recalibrated_bams, recalibrated_bams_bcftools, recalibrated_bams_qc
  
  script:
  // TODO: add BED file?
  """
  /gatk/gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -bqsr ${recalibration_table} \
      -I ${bam} \
      -O ${name}.recal.bam \
      -R ${fasta}
  """
}

/*--------------------------------------------------
  Run BAM QC using Qualimap after recalibration
---------------------------------------------------*/

process runBamQCrecalibrated {
    tag "$bam"
    container 'maxulysse/sarek:latest'

    input:
    set val(name), file(bam), file(index) from recalibrated_bams_qc

    output:
    file("${name}_recalibrated") into bamqc_recalibrated_report

    when: !params.skip_multiqc

    script:
    """
    qualimap \
      bamqc \
      -bam ${bam} \
      --paint-chromosome-limits \
      --genome-gc-distr HUMAN \
      -nt ${task.cpus} \
      -skip-duplicated \
      --skip-dup-mode 0 \
      -outdir ${name}_recalibrated \
      -outformat HTML
    """
}

haplotypecaller = recalibrated_bams.combine(haplotypecaller_ref)

/*--------------------------------------------------
  HaplotypeCaller per-sample
---------------------------------------------------*/

process haplotypeCaller {
  tag "${name}"
  publishDir "${params.outdir}/haplotypeCaller", mode: 'copy'
  
  container 'broadinstitute/gatk:4.0.4.0'
  memory "14.GB"
  cpus 4

  input:
  set val(name), file(bam), file(bai), file(fasta), file(fai), file(dict) from haplotypecaller

  output:
  set val(name), file("${name}.GATK.vcf"), file("${name}.GATK.vcf.idx") into haplotypecaller_vcf

  script:
  // TODO: add BED, intervals & dbSNP?
  """
  /gatk/gatk --java-options "-Xmx4g" \
      HaplotypeCaller \
      -R ${fasta} \
      -I ${bam} \
      -O ${name}.GATK.vcf 
  """
}

bcftools = recalibrated_bams_bcftools.combine(bcftools_ref)

/*--------------------------------------------------
  bcftools Mpileup variant calling
---------------------------------------------------*/

process bcftoolsMpileup {
  tag "${name}"
  publishDir "${params.outdir}/bcftoolsMpileup", mode: 'copy'
  
  container "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
  memory "14.GB"
  cpus 4

  input:
  set val(name), file(bam), file(bai), file(fasta), file(fai), file(dict) from bcftools

  output:
  set val(name), file("${name}.SAM.vcf") into bcftools_vcf

  script:
  """
  bcftools mpileup \
      --max-depth 10000 \
      --max-idepth 10000 \
      --annotate "FORMAT/AD,FORMAT/DP" \
      --fasta-ref ${fasta} \
      --ignore-RG \
      --no-BAQ \
      ${bam} | bcftools call -Ov -mv \
          -o ${name}.SAM.vcf
  """
}

vcfs = haplotypecaller_vcf.combine(bcftools_vcf, by: 0)
annovar_consensus = vcfs.combine(annotation)

/*--------------------------------------------------
  annotate with annovar
---------------------------------------------------*/

process annovarConsensus {
  tag "${name}"
  publishDir "${params.outdir}/annotation", mode: 'copy'

  container "bioinfochrustrasbourg/annovar:2018Apr16"
  memory "4.GB"

  input:
  set val(name), file(gatk_vcf), file(gatk_vcf_index), file(sam_vcf), file(annotation) from annovar_consensus

  output:
  set file("${name}/${name}.GATK.${params.ref_name}_multianno.vcf"), file("${name}/${name}.GATK.${params.ref_name}_multianno.txt"),
    file("${name}/${name}.SAM.${params.ref_name}_multianno.vcf"), file("${name}/${name}.SAM.${params.ref_name}_multianno.txt") into annotated_results

  script:
  """
  mkdir ${name}

  table_annovar.pl ${gatk_vcf} ${annotation} \
    -buildver ${params.ref_name} \
    -outfile ${name}/${name}.GATK \
    -remove \
    -protocol ${params.annovar_protocols} \
    -operation ${params.annovar_operation} \
    -nastring . -vcfinput

  table_annovar.pl ${sam_vcf} ${annotation} \
    -buildver ${params.ref_name} \
    -outfile ${name}/${name}.SAM \
    -remove \
    -protocol ${params.annovar_protocols} \
    -operation ${params.annovar_operation} \
    -nastring . -vcfinput
  """
}

/*--------------------------------------------------
  bcftools stats to get summary statistics for VCFs
---------------------------------------------------*/

process bcftoolsStats {
  tag "${name}"
  
  container "quay.io/biocontainers/bcftools:1.9--h4da6232_0"

  input:
  set file(gatk_vcf), file(gatk), file(sam_vcf), file(sam) from annotated_results

  output:
  set file("bcfstats_${gatk_vcf}.txt"), file("bcfstats_${sam_vcf}.txt") into bcftools_stats

  when: !params.skip_multiqc

  script:
  """
  bcftools stats $gatk_vcf > bcfstats_${gatk_vcf}.txt
  bcftools stats $sam_vcf > bcfstats_${sam_vcf}.txt
  """
}

/*--------------------------------------------------
  Make output report using MultiQC
---------------------------------------------------*/

process multiqc {
  tag "multiqc_report.html"

  publishDir "${params.outdir}/MultiQC", mode: 'copy'
  container 'ewels/multiqc:v1.7'

  input:
  file fastqc from fastqc_results.collect()
  file bamqc_mapped from bamqc_mapped_report.collect()
  file bamqc_recalibrated from bamqc_recalibrated_report.collect()
  file bcftools_stats from bcftools_stats.collect()

  output:
  file("*") into viz

  when: !params.skip_multiqc

  script:
  """
  multiqc . -m fastqc -m qualimap -m bcftools 
  """
}