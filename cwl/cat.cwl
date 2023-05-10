#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  Merge a set of files into file using the cat utility.

requirements:
 - class: InlineJavascriptRequirement

baseCommand: cat

inputs:
  output_file:
    type: string
    doc: |
      Name of final file.

  input_files:
    type: File[]
    inputBinding:
      position: 4
    doc: |
      Array of files to merge. 

stdout: $(inputs.output_file)

outputs:
  combined_file:
    type: File
    outputBinding:
      glob: $(inputs.output_file)