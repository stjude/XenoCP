#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  Merge a set of sorted bam files into one bam file, 
  mark duplicates, and index. Also calculates flagstate
  and md5. 

requirements:
 - class: InlineJavascriptRequirement

baseCommand: merge_markdup_index.sh

inputs:
  skip_dup:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --no-markdup
    doc: |
      Flag as true in order to skip
      the mark duplicates step.

  output_bam:
    type: string
    inputBinding:
      position: 2
    doc: |
      Name of final output bam file.
      Must have '.bam' extension. 

  input_bams:
    type: File[]
    inputBinding:
      position: 3
    doc: |
      Array of bam files to merge. 

outputs:
  final_bam:
    type: File
    secondaryFiles: 
      - .bai
      - .md5
    outputBinding:
      glob: $(inputs.output_bam)
  
  flagstat:
    type: File
    outputBinding:
      glob: |
        ${ 
          return inputs.output_bam.replace(".bam", ".flagstat.txt");
        }

#  markdup:
#    type: File
#    outputBinding:
#      glob: |
#        ${
#          var i = inputs.output_bam.split('.').slice(0, -1).join('.');
#          if (inputs.skip_dup) {
#            return [];
#          } else {
#            return i + ".markdup.txt";
#          }
#        }


