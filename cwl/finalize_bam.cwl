#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
doc: |
  Clean duplicate BAM header records from STAR aligned BAMs.
  Add @PG record for XenoCP Also calculates flagstat
  and md5. 

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      finalize_bam:
        specs: ["finalize_bam.sh"]

baseCommand: finalize_bam.sh

inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
    doc: |
      BAM file to clean. 
  output_bam:
    type: string
    inputBinding:
      position: 2
    doc: |
      Name of final output bam file.
      Must have '.bam' extension. 
  ref_db_prefix: 
    type: string
    inputBinding:
      position: 3
  aligner:
    type:
      type: enum
      symbols: ["bwa aln", "bwa mem", "star"]
      name: aligner
    inputBinding:
      position: 4
  original_bam_name:
    type: string
    inputBinding:
      position: 5

outputs:
  final_bam:
    type: File
    secondaryFiles: 
      - .md5
      - .bai
    outputBinding:
      glob: $(inputs.output_bam)