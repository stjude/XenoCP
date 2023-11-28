#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: star_onlymapped.sh

requirements:
  ResourceRequirement:
    ramMin: 30000
    coresMin: 2

hints:
  SoftwareRequirement:
    packages:
      star:
        specs: ["STAR"]
      tweak_sam:
        specs: ["java.sh org.stjude.compbio.sam.TweakSam"]
  DockerRequirement:
    dockerPull: "ghcr.io/stjude/xenocp:latest"

inputs:
  ref_db_prefix: 
    type: string
    inputBinding:
        position: 4
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
    inputBinding:
        position: 1

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
  star_onlymapped.sh
    Performs STAR to align to a reference, converts to bam, and
    filters out all records for unmapped reads.

    $1 = reference db
    $2 = input fastq file
    $3 = output bam

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
