#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  bam: 
    type: File
    label: Aligned sequences in BAM format
    secondaryFiles:
     - .bai

outputs:
  split_bams:
    type: File[]
    outputSource: [by_chrom/out_bam, mismatch/out_bam]
    linkMerge: merge_flattened
  unmapped_bams:
    type: File[]
    outputSource: [unmapped/out_bam]

steps:
  # Step01: extract chromosome information from input bam
  get_chroms: 
    in: 
      bam: bam
    out: [chroms]
    run: get_chroms.cwl  

  # Step 02a: extract reads with mates mapped to other chromosomes
  mismatch: 
    in:
      bam: bam
    out: [out_bam]
    run: 
      class: CommandLineTool
      stdout: other.bam
      inputs: 
        bam: 
          type: File
          inputBinding: 
            position: 0
            prefix: -i 
      arguments: ["-x", "-O", "queryname", "-o", "other.bam"]
      outputs: 
        out_bam: 
          type: File
          outputBinding: 
            glob: "*.bam"
      baseCommand: [java.sh, org.stjude.compbio.sam.TweakSam] 

  # Step02b: extract reads mapped to each chromosome
  by_chrom:
    in:
      chroms:  
        source: get_chroms/chroms
      bam: bam
    out: [out_bam]
    scatter: chroms
    run:
      class: CommandLineTool
      inputs:
        chroms:
          type: string
          inputBinding:
            position: 1
            prefix: -c
        bam:
          type: File
          inputBinding:
            position: 0
            prefix: -i
      arguments: ["-X", "-O", "queryname", "-o", "$(inputs.chroms).bam"]
      outputs: 
        out_bam:
         type: File
         outputBinding: 
           glob: "$(inputs.chroms).bam"
      baseCommand: [java.sh, org.stjude.compbio.sam.TweakSam]

  # Step 02c: extract unmapped reads
  unmapped:
    in:
      bam: bam
    out: [out_bam]
    run:
      class: CommandLineTool
      stdout: unmapped.bam
      inputs:
        bam:
          type: File
          inputBinding:
            position: 0
            prefix: -i
      arguments: ["-n", "-O", "coordinate", "-o", "unmapped.bam"]
      outputs:
        out_bam:
          type: File[]
          outputBinding:
            glob: "*.bam"
      baseCommand: [java.sh, org.stjude.compbio.sam.TweakSam]

