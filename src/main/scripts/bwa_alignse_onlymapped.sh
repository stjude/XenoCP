#!/bin/bash
# Performs bwa aln and bwa samse to align to a reference, converts to bam, and
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

set -e -x

# bwa aln
cmd="bwa aln $PREFIX $FASTQ"
sai=$SCRATCH/`basename $BAM .bam`.sai
echo $cmd \> $sai
$cmd > $sai

# bwa samse | TweakSam
cmd="bwa samse $PREFIX $sai $FASTQ | java.sh org.stjude.compbio.sam.TweakSam -V SILENT -G 4 -O queryname -o $BAM"
echo "$cmd"
eval $cmd

# Remove scratch dir
rm -rf $SCRATCH

# Write end date
date
