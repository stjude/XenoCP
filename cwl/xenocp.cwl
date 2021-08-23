#!/usr/bin/env cwl-runner
cwlVersion: v1.2
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
  aligner:
    type:
      type: enum
      symbols: ["bwa aln", "bwa mem", "star"]
      name: aligner
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
    type: File
    outputSource: merge_contam_list/combined_file
  output_tie_bam:
    type: File
    outputSource: merge_tie_bam/final_bam

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

  # Step03a: map extracted reads to the contamination genome with bwa aln
  mapping-bwa-aln:
    run: bwa_alignse_onlymapped.cwl
    when: $(inputs.aligner == "bwa aln")
    in:
      aligner: aligner
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

  # Step03b: map extracted reads to the contamination genome with bwa mem
  mapping-bwa-mem:
    run: bwa_mem_onlymapped.cwl
    when: $(inputs.aligner == "bwa mem")
    in:
      aligner: aligner
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

  # Step03c: map extracted reads to the contamination genome with STAR
  mapping-star:
    run: star_onlymapped.cwl
    when: $(inputs.aligner == "star")
    in:
      aligner: aligner
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
      contam_bams:
        source: [mapping-bwa-aln/bam, mapping-bwa-mem/bam, mapping-star/bam]
        linkMerge: merge_flattened
        pickValue: all_non_null
    scatter: [input_bam, contam_bams]
    scatterMethod: dotproduct
    out: [contam_list, output_tie_bam]

  # Step05a: clean the original bam by setting the contamination reads to be unmapped
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

  # Step05b: sort tie BAMs prior to merge
  sort-bams:
    run: bio-cwl-tools:picard/picard_SortSam.cwl
    in:
      alignments: contamination/output_tie_bam
      sort_order:
        valueFrom: $("coordinate")
      validation_stringency:
        valueFrom: $("SILENT")
    scatter: [alignments]
    scatterMethod: dotproduct
    out: [sorted_alignments]

  # Step06a: merge split bams, index and mark duplicates
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

  # Step06b: merge tie bams and index
  merge_tie_bam:
    run: merge_markdup_index.cwl
    in:
      input_bams:
        source: [cleanse/cleaned_bam, split/unmapped]
        linkMerge: merge_flattened
      output_bam:
        source: bam
        valueFrom: ${return self.nameroot + ".xenocp.bam"}
      n_threads: n_threads
      skip_dup:
        valueFrom: $(true)
    out: [final_bam, flagstat]

  # Step07: Combine contam lists
  merge_contam_list:
    run: cat.cwl
    in:
      input_files: contamination/contam_list
      output_file:
        source: bam
        valueFrom: ${return self.nameroot + ".contam.txt"}
    out: [combined_file]

  # Step08: QC the merged bam
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
  bio-cwl-tools: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
