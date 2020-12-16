#!/bin/bash
# Performs bwa mem to align to a reference, converts to bam, and
# filters out all records for unmapped reads.
#
# $1 = reference db prefix
# $2 = input fastq file
# $3 = output bam

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PREFIX=$1
FASTQ=$2
BAM=$3

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $BAM`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Create temp dir to hold intermediate files
SCRATCH=`mktemp -d` || ( echo "Could not get scratch dir" >&2 ; exit 1 )

# bwa mem
cmd="bwa mem $PREFIX $FASTQ | java.sh org.stjude.compbio.sam.TweakSam -V SILENT -G 4 -o $BAM"

echo $cmd 
set -e -x
eval $cmd

# Remove scratch dir
rm -rf $SCRATCH

# Write end date
date
