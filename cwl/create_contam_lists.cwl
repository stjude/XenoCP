#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [java.sh, org.stjude.compbio.xenocp.CreateContamLists]

hints:
  SoftwareRequirement:
    packages:
      create_contam_list:
        specs: [ "java.sh org.stjude.compbio.xenocp.CreateContamLists" ]
  DockerRequirement:
    dockerPull: "ghcr.io/stjude/xenocp:latest"

inputs:
  input_bam:
    type: File
    label: Aligned sequences in SAM/BAM format
    inputBinding:
      position: 1
      prefix: -i
  output_contam_list:
    type: string
    label: Must be an output contam list file name, not an absolute path
    inputBinding:
      position: 2
      prefix: -o
  tie_bam:
    type: string
    label: Must be an output tie bam file name, not an absolute path
    inputBinding:
      position: 3
      prefix: -t
  stringency:
    type: string?
    default: SILENT
    inputBinding:
     position: 4
     prefix: -V
  contam_bams:
    type: File
    label: comma separated list of bams mapped to contamination bams
    inputBinding:
      position: 5

outputs:
  contam_list:
    type: File
    outputBinding:
      glob: $(inputs.output_contam_list)
  output_tie_bam:
    type: File
    outputBinding:
      glob: $(inputs.tie_bam)

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
  java.sh org.stjude.compbio.xenocp.CreateContamLists
    usage: java org.stjude.compbio.xenocp.CreateContamLists [OPTION] CONTAM_BAM1,...
                Looks for contamination in a sam/bam file by finding a
                better mapping in any of a list of contamination bams.
                The main output is a list of read names that are called as
                contamination.  You may also use -t to write records that
                were called as 'ties' to an output sam/bam file for
                review.
                This will operate on single- or paired-end data, but the
                input and contam bams must be 'logically' both single or
                or paired.  A read is considered 'logically' paired if it
                is actually paired or if it is single end but has a
                suffix of [non-alphanumeric char][12] that indicates
                read number in the pair.
                CAVEAT: All bams must have the NM tag, containing edit
                distance!
    -i <arg>   input file to use if not stdin
    -o <arg>   output contam reads file if not stdout
    -t <arg>   output sam/bam to which to write ties
    -V <arg>   validation stringency: STRICT, LENIENT, or SILENT (default: SILENT)

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
