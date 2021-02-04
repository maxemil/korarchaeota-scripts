#!/usr/bin/env nextflow

params.fasta = ""
params.outdir = ""

fasta_input = Channel.fromPath(params.fasta)

process runEMapper {
  input:
  file fasta from fasta_input

  output:
  file "${fasta.simpleName}.emapper.seed_orthologs" into seed_orthologs
  file "${fasta.simpleName}.emapper.annotations" into annotations

  publishDir "${params.outdir}", mode: 'copy'
  tag {"${fasta.baseName}"}
  cpus 2

  script:
  """
  # emapper2
  /local/two/Software/eggnog-mapper2/emapper.py -o ${fasta.simpleName} \
                                                -i $fasta \
                                                -m hmmer \
                                                -d archaea \
                                                --cpu ${task.cpus} \
                                                --tax_scope 2157 
  """
}
