#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: bwa_mem_onlymapped.sh

hints:
  SoftwareRequirement:
    packages:
      bwa:
        specs: ["bwa mem"]
      tweak_sam:
        specs: ["java.sh org.stjude.compbio.sam.TweakSam"]
  DockerRequirement:
    dockerPull: "ghcr.io/stjude/xenocp:latest"

requirements:
  ResourceRequirement:
    ramMin: 15000
    coresMin: 1

inputs:
  ref_db_prefix: 
    type: string
    inputBinding:
        position: 1
        valueFrom: |
          ${
            return inputs.index.path + "/" + self;
          }
  input_fastq:
    type: File
    inputBinding:
        position: 2
  output_bam:
    type: string
    label: Must be an output bam file name, not an absolute path
    inputBinding:
        position: 3
  index:
    type: Directory

outputs:
    bam:
      type: File
      outputBinding:
          glob: $(inputs.output_bam)

s:author:
  class: s:Person
  s:name: Andrew Thrasher
  s:email: mailto:Andrew.Thrasher@stjude.org
  s:worksFor:
  - class: s:Organization
    s:name: St.Jude Children's Research Hospital
    s:location: 262 Danny Thomas Place Memphis, TN 38105
    s:department:
    - class: s:Organization
      s:name: Computational Biology

doc: |
  bwa_mem_onlymapped.sh
    Performs bwa mem to align to a reference, converts to bam, and
    filters out all records for unmapped reads.

    $1 = reference db prefix
    $2 = input fastq file
    $3 = output bam

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
