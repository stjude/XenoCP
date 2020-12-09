#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  bam: 
    type: File
    label: Aligned sequences in BAM format
    secondaryFiles:
     - .bai
  ref_db_prefix:
    type: string
    label: contamination genome reference db prefix
  # See doc in split_sam.cwl for the meaning of the following arguments
  suffix_length:
    type: int?
    default: 4
  keep_mates_together:
    type: boolean?
    default: True
  validation_stringency:
    type: string?
    default: SILENT
  output_prefix:
    type: string?
    default: xenocp-
  output_extension:
    type: string?
    default: bam
  n_threads: 
    type: int?
    default: 1

outputs:
  output_bam:
    type: File
    outputSource: finish/final_bam
  flagstat:
    type: File
    outputSource: finish/flagstat
  contam_list:
    type:
      type: array
      items: File
    outputSource: contamination/contam_list
  output_tie_bam:
    type:
      type: array
      items: File
    outputSource: contamination/output_tie_bam

steps:
  # Step01: extract chromosome information from input bam
  split:
    run: extract.cwl
    in: 
      bam: bam
    out: [split_bams, unmapped]
  
  # Step02: extract mapped reads and convert to fastq
  mapped-fastq:
    run: view_awk_picard.cwl
    in:
      input_bam: 
        source: split/split_bams
        linkMerge: merge_flattened 
      output_fastq:
        valueFrom: $(inputs.input_bam.nameroot).fastq   # here "inputs" refers to inputs in view_awk_picard.cwl
    scatter: [input_bam]
    out: [fastq]

  # Step03: mapped extacted reads to the contamination genome
  mapping:
    run: bwa_alignse_onlymapped.cwl
    in:
      ref_db_prefix: ref_db_prefix
      input_fastq: mapped-fastq/fastq
      output_bam:
        valueFrom: $(inputs.input_fastq.nameroot).contam.bam
    scatter: [input_fastq]
    out: [bam]
    hints:
      ResourceRequirement:
        ramMin: 4800
        coresMin: 1

  # Step04: find contamination reads
  contamination:
    run: create_contam_lists.cwl
    in:
      input_bam:
        source: split/split_bams 
      output_contam_list: 
        valueFrom: $(inputs.input_bam.nameroot).contam.txt
      tie_bam:
        valueFrom: $(inputs.input_bam.nameroot).tie.bam
      contam_bams: mapping/bam
    scatter: [input_bam, contam_bams]
    scatterMethod: dotproduct
    out: [contam_list, output_tie_bam]
  
  # Step05: clean the original bam by setting the contamination reads to be unmapped
  cleanse:
    run: tweak_sam.cwl
    in:
      input_bam: 
        source: split/split_bams 
      output_bam:
        valueFrom: $(inputs.input_bam.nameroot).cleaned.bam
      unmap_reads: contamination/contam_list
    scatter: [input_bam, unmap_reads]
    scatterMethod: dotproduct
    out: [cleaned_bam]

  # Step06: merge split bams, index and mark duplicates
  finish:
    run: merge_markdup_index.cwl
    in:
      input_bams: 
        source: [cleanse/cleaned_bam, split/unmapped]
        linkMerge: merge_flattened
      output_bam:
        source: bam
        valueFrom: ${return self.nameroot + ".xenocp.bam"}
      n_threads: n_threads
    out: [final_bam, flagstat]
  
  # Step07: QC the merged bam
  finalqc:
    run: qc_bam.cwl
    in:
      bam: finish/final_bam
      flagstat: finish/flagstat
    out: []

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
    Xenograft sample cleaning pipeline

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
