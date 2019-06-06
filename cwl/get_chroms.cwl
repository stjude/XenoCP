#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

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
      #glob: $(stdout)
      loadContents: true
      outputEval: $(self[0].contents.split("\n").slice(0,-1))
    type: string[]
