#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [java.sh, org.stjude.compbio.sam.SplitSam]

hints:
  SoftwareRequirement:
    packages:
      SplitSam:
        specs: [ "SplitSam.java" ]

inputs:
  suffix_length:
    type: int?
    inputBinding:
      position: 2
      prefix: -a
  add_rgid:
    type: boolean?
    inputBinding:
      position: 6
      prefix: --add-unique-rgid-to-records
  num_bucket:
    type: int?
    inputBinding:
      position: 4
      prefix: -b
  split_by_ref_name:
    type: boolean?
    inputBinding:
      position: 7
      prefix: -c
  input_sam_file:
    type: File?
    inputBinding:
      position: 1
      prefix: -i
  seq_len_divisor:
    type: string?
    inputBinding:
      position: 8
  keep_mates_together:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -m
  strips_mate_suffixes:
    type: boolean?
    inputBinding:
      position: 9
      prefix: -M
  split_by_arg:
    type: string?
    inputBinding:
      position: 10
      prefix: -n
  filter_out_non_primary:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -p
  split_by_read_group:
    type: boolean?
    inputBinding:
      position: 12
      prefix: -r
  validation_stringency:
    type: string?
    inputBinding:
      position: 5
      prefix: -V
  output_prefix:
    type: string
    inputBinding:
      position: 13
  output_extension:
    type: string
    inputBinding:
      position: 14

outputs:
  split_bams:
    type:
      type: array
      items: File
    outputBinding:
      glob: $(inputs.output_prefix)*.bam

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
  CompBio in-house program to split a sam/bam file into smaller pieces
    Usage: java org.stjude.compbio.sam.SplitSam [OPTION]... PREFIX EXT
            Splits a sam/bam file into smaller pieces. Output is
            named with prefix PREFIX and extension EXT, with
            distinguishing parts in between. The extension
            determines the output format. PREFIX may contain path.
            You may split by reference name, read group, and/or by
            number of reads.
            If your data is paired-end, and you would like to preserve
            sort order while keeping mates in the same file, then you
            can use -m with -b (and optionally -l). In this mode, all
            inputs are divided by read name into B buckets. You can
            specify the number of buckets using -b, or compute by
            taking chromosome length divided by the -l value. If you use -l,
            you still need -b for the no-ref case.
    -a <arg>                          suffix length for -n default 3
    --add-unique-rgid-to-records      add RGID to records where missing, if
                                      and only if there is exactly one RG in header
    -b <arg>                          number of buckets for -m splitting
    -c                                split by ref name (usually chromosome)
    -i <arg>                          input sam/bam file if not stdin
    -l <arg>                          seq len divisor to compute b on the fly
    -m                                keep mates together
    -M                                strips mate suffixes /1 and /2
    -n <arg>                          split every 'arg' records
    -p                                filter out non-primary records
    -r                                split by read group
    -V <arg>                          validation stringency: STRICT, LENIENT,
                                      or SILENT (default: SILENT)

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
