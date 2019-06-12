#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      bam_to_chr:
        specs: ["bam_to_chrs.sh"]


baseCommand: bam_to_chrs.sh

requirements: 
  InlineJavascriptRequirement: {} 

stdout: $(inputs.bam.nameroot).chrs

inputs:
  bam: 
    type: File
    inputBinding: 
      position: 1
outputs:
  chroms: 
    outputBinding: 
      glob: "*.chrs"
      loadContents: true
      outputEval: $(self[0].contents.split("\n").slice(0,-1))
    type: string[]
