#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: sort_flagstat.sh

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      sort_flagstat:
        specs: ["sort_flagstat.sh"]
      picard:
        specs: ["java.sh picard.cmdline.PicardCommandLine SortSam"]
        version: ["2.6.0"]
      samtools:
        specs: ["samtools flagstat"]
        version: ["1.3.1"]

inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
    doc: Aligned sequences in BAM format
  sort_order: 
    type: string
    inputBinding:
      position: 2
  output_bam:
    type: string
    inputBinding:
      position: 3
    doc: Must be an output bam file name, not an absolute path
  other_params:
    type: string?
    inputBinding:
      position: 4

outputs:
  sort_bam:
    type: File
    outputBinding:
      glob: $(inputs.output_bam)
  flagstat:
    type: File
    outputBinding:
      glob: $(inputs.output_bam.split('.')[0]).flagstat.txt

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
  Performs Picard SortSam and then flagstat
    Flagstat file will have same name as output BAM, but with .flagstat.txt in
    place of .bam
      Parameters:
      $1 = input file
      $2 = sort order
      $3 = output BAM file
      $4,... = other parameters are sent to SortSam

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
