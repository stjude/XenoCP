#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: view_awk_picard.sh

hints:
  SoftwareRequirement:
    packages:
      sort_flagstat:
        specs: ["sam_to_single.awk"]
      picard:
        specs: ["picard.cmdline.PicardCommandLine SamToFastq"]
        version: ["2.6.0"]
      samtools:
        specs: ["samtools view"]
        version: ["1.3.1"]

inputs:
  input_bam:
    type: File
    label: Aligned sequences in BAM format
    inputBinding:
      position: 1
    streamable: true
  output_fastq:
    type: string
    label: Must be an output fastq file name, not an absolute path
    inputBinding:
      position: 2
    streamable: true

outputs:
  fastq:
    type: File
    outputBinding:
      glob: $(inputs.output_fastq)

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
  view_awk_picard.sh
    A wapper of samtools view, sam_to_single.awk and picard SamToFastq for XenoCP pipeline to
    extract mapped reads and convert to fastq.
    Parameters:
    $1 = input aligned bam file
    $2 = output fastq file

$namespaces:
  s: http://schema.org

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
