#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [java.sh, org.stjude.compbio.sam.TweakSam]

hints:
  SoftwareRequirement:
    packages:
      create_contam_list:
        specs: [ "java.sh org.stjude.compbio.sam.TweakSam" ]

inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
      prefix: -i
  output_bam:
    type: string
    label: Must be an output cleaned bam file name, not an absolute path
    inputBinding:
      position: 2
      prefix: -o
  unmap_reads:
    type: File
    inputBinding:
      position: 3
      prefix: -u
  sort_order:
    type: string?
    default: coordinate
    inputBinding:
      position: 4
      prefix: -O
  stringency:
    type: string?
    default: SILENT
    inputBinding:
     position: 5
     prefix: -V

outputs:
  cleaned_bam:
    type: File
    outputBinding:
      glob: $(inputs.output_bam)

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
  java.sh org.stjude.compbio.sam.TweakSam
    usage: java -jar TweakSam.jar [OPTION]...
                Copies a sam/bam file while performing any/all of the following tweaks among others:
                1. Format conversion.  Output format is determined by extension.  Stdout output is always SAM.
                2. Sorting for almost sorted data.  If the data is almost sorted, then sorting may be fixed using the -f option.
                  A buffer size may be specified with -b.  If any record is further away in the input than (buffer size), then
                  the command will fail.  The sort order to use is given by -O.
                3. Sort order declaration.  If the data is already sorted but the order is not declared, or if you are fixing
                  sorting using -f, you may use -O to declare the correct sort order.  This will be set in the output header.
                4. Extraction by sequence name.  Specify a reference sequence name using -c or specify -n to extract records
                  with no reference sequence name.  These are only valid for coordinate-sorted, indexed bam files5. CIGAR overflow
                  detection.  Detects overflow of the numeric part of a CIGAR element and soft clips the alignment at that position.
                6. Random sampling.  Specify a probability using -r to randomly sample from the input bam.  Mate pairs should
                  NOT be broken up, but sorting is preserved.
                7. Alignment and allele selection.  Specify a position using -c and -p to extract reads aligned at that
                  position.  Use -a to specify the allele at that position, or -L and/or -R to specify an interval over
                  which there must be contiguous alignment.  Use -M to require a certain match rate on either side of the -p position.
                8. Flag filtering using -g and -G like samtools view -f and -F
                9. Others (see arg list)
    -5                          create MD5 file while writing output
    -a <arg>                    show only records with the given allele at the position specified with -p 
                                (valid values: C, G, T, A; requires -p
    -C                          exclude consecutive duplicates in queryname-sorted input (these are not marked dups, 
                                they are dups in read name/pairing flags)
    -c <arg>                    only include records with the given reference name (requires indexed bam input)
    -D                          exclude marked duplicates
    -d                          clear duplicate flag
    --filter-read-id-dups       filter out records that are duplicates based on read id (read name + pairing
                                flags).  When there are multiple records for the same read, only the first is
                                written. Note: significantly increases memory requirements.
    -G <arg>                    show only records with the specified flag(s) unset
    -g <arg>                    show only records with the specified flag(s) set
    -h <arg>                    replace header with one from the specified file
    -i <arg>                    input sam/bam file if not stdin
    -L <arg>                    requires -p; indicates that the read must also align for
                                n bases to the left of p, for a total of n+1 aligned bases
    -l <arg>                    only include records with read length in range; formats are N >N <N [M,N] (M,N) [M,N) (M,N]
    -M <arg>                    minimum match rate in the alignment block on either side of the position indicated by -p.
                                The block is split at the -p position, and each part of the block is tested
                                separately; both parts must pass for the read to be included.
    -m <arg>                    maximum number of reads to return
    -n                          only include records with no reference name (requires indexed bam input)
    --normalize-pairs           Aggressively works to normalize read pairs, regardless of input sorting.  May
                                dramatically increase memory requirements.
    -O <arg>                    sort order: coordinate, queryname, or unsorted
    -o <arg>                    output sam/bam file if not stdout
    -p <arg>                    show only records that have a base aligned to the given position (requires -c)
    -q                          filter out records marked as 'duplicate' or 'fails vendor quality checks'
    -R <arg>                    requires -p; indicates that the read must also align for n bases to the right of p, 
                                for a total of n+1 aligned bases
    -r <arg>                    randomly sample records with the given probability of including any one record. Range is (0.0,1.0).
    --set-sort-order <ORDER>    output sort order to declare: queryname, coordinate, or unsorted.  Input must
                                already be sorted in the proper order; if it is not, then use -O
    -U                          strip all alignment info; make it a UBAM
    -u <arg>                    set reads named in the file to unmapped; add XU tag
    -V <arg>                    validation stringency: STRICT, LENIENT, or SILENT (default: SILENT)
    -X                          exclude paired records whose mates are mapped to a different chromosome
    -x                          only include paired records whose mates are mapped to a different chromosome
    -y <arg>                    assign read groups according to mapping in file; the file is tab-delimited text with a
                                regex and an arbitrary identifier.  If the read name matches the regex, then the read
                                group ID is set to the ID indicated in the file

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
