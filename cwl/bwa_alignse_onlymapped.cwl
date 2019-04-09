#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: bwa_alignse_onlymapped.sh

hints:
  SoftwareRequirement:
    packages:
      bwa:
        specs: ["bwa aln", "bwa samse"]
      tweak_sam:
        specs: ["java.sh org.stjude.compbio.sam.TweakSam"]

inputs:
  ref_db_prefix: 
    type: string
    inputBinding:
        position: 1
  input_fastq:
    type: File
    inputBinding:
        position: 2
  output_bam:
    type: string
    label: Must be an output bam file name, not an absolute path
    inputBinding:
        position: 3

outputs:
    bam:
      type: File
      outputBinding:
          glob: $(inputs.output_bam)

s:author:
  class: s:Person
  s:name: Liang Ding
  s:email: mailto:Liang.Ding@stjude.org
  s:worksFor:
  - class: s:Organization
    s:name: St.Jude Children's Research Hospital
    s:location: 262 Danny Thomas Place Memphis, TN 38105
    s:department:
    - class: s:Organization
      s:name: Computational Biology

doc: |
  bwa_alignse_onlymapped.sh
    Performs bwa aln and bwa samse to align to a reference, converts to bam, and
    filters out all records for unmapped reads.

    $1 = reference db prefix
    $2 = input fastq file
    $3 = output bam

$namespaces:
  s: http://schema.org

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
