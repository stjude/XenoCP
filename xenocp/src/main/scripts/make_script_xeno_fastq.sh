#!/bin/bash
# Writes a script to extract fastq files from a set of bam files for use with
# the Xenograft purification pipeline.
#
# Only bams not called noref* will be included, and only mapped reads are
# extracted.  Paired end reads are "cast" to single-end by adding a .M suffix
# to the read names (where M is 1 or 2).
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = input directory containing bams
# $2 = output directory for fastqs

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
IN_DIR=$1
OUT_DIR=$2

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# One job per bam
for bam in $IN_DIR/*.bam
do
  bn=`basename $bam .bam`
  echo samtools view -h -F 4 $bam \| sam_to_single.awk -v delim=. \| java.sh picard.cmdline.PicardCommandLine SamToFastq INPUT=/dev/stdin FASTQ=$OUT_DIR/$bn.fq VALIDATION_STRINGENCY=SILENT
done
