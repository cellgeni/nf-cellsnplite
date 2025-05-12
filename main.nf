#!/usr/bin/env nextflow

process MODE_2B {

  publishDir "${launchDir}/cellsnplite-results-${params.project_tag}/step1", mode: 'copy'

  input:
  tuple val(sampleID), path(bam), path(barcodes), path(bai)

  output:
  tuple val(sampleID), path(bam), path(barcodes), path(bai), path("${sampleID}/*.vcf.gz*")
  path("${sampleID}/*.tsv")
  path("${sampleID}/*.mtx")

  script:
  """
  cellsnp-lite -s ${bam} -O ${sampleID} -p 10 --minMAF ${params.minMAF} --minCOUNT ${params.minCOUNT} --cellTAG None --UMItag ${params.UMItag} ${params.refseq ? "-f ${params.refseq_path}" : ""} --gzip
  tabix -p vcf ${sampleID}/cellSNP.base.vcf.gz 
  """
}

process MODE_1A {

  publishDir "${launchDir}/cellsnplite-results-${params.project_tag}/step2", mode: 'copy'

  input:
  tuple val(sampleID), path(bam), path(barcodes), path(bai), path(vcf), path(vcf_tbi)

  output:
  path("${sampleID}/*.vcf.gz")
  path("${sampleID}/*.tsv")
  path("${sampleID}/*.mtx")

  script:
  """
  cellsnp-lite -s ${bam} -b ${barcodes} -O ${sampleID} -R ${vcf} -p 10 --minMAF ${params.minMAF} --minCOUNT ${params.minCOUNT} --cellTAG ${params.cellTAG} --UMItag ${params.UMItag} ${params.genotype ? '--genotype' : ''} --gzip
  """
}

process MODE_2A {

  publishDir "${launchDir}/cellsnplite-results-${params.project_tag}/", mode: 'copy'

  input:
  val(sampleID)
  path(bam)
  path(barcodes)

  output:
  path("${sampleID}/*.vcf.gz")
  path("${sampleID}/*.tsv")
  path("${sampleID}/*.mtx")

  script:
  """
  cellsnp-lite -s ${bam} -b ${barcodes} -O ${sampleID} -p 10 --minMAF ${params.minMAF} --minCOUNT ${params.minCOUNT} ${params.genotype ? '--genotype' : ''} --gzip
  """
}



workflow '2step' {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv( header: ['sampleID', 'bam', 'barcodes'],  sep: '\t')
       .map { row -> tuple(row.sampleID, file(row.bam), file(row.barcodes), file(row.bam + '.bai')) }
       .set{ row }
  MODE_2B(row)
  MODE_1A(MODE_2B.out[0])
}

workflow '1step' {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv( header: ['sampleID', 'bam', 'barcodes'],  sep: '\t')
       .map { row -> tuple(row.sampleID, file(row.bam), file(row.barcodes), file(row.bam + '.bai')) }
       .set{ row }
  MODE_2A(row)
}

workflow '1a-only' {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv(header: ['sampleID', 'bam', 'barcodes', 'vcf'], sep: '\t')
       .map { row -> tuple(row.sampleID, file(row.bam), file(row.barcodes), file(row.bam + '.bai'), file(row.vcf)) }
       .set { row }
  MODE_1A(row)
}

workflow '2b-only' {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv( header: ['sampleID', 'bam', 'barcodes'],  sep: '\t')
       .map { row -> tuple(row.sampleID, file(row.bam), file(row.barcodes), file(row.bam + '.bai')) }
       .set{ row }
  MODE_2B(row)
}