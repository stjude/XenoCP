#!/bin/bash
# Writes a script to run bwa_alignse_onlymapped.sh on a set of fastq files in a
# directory
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = genome
# $2 = input directory containing fastqs
# $3 = output directory for bams

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOME=$1
IN_DIR=$2
OUT_DIR=$3

# Import configuration variable for whole genome bwa database
. import_config.sh genome $GENOME WG_BWA_DB

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# One job per fastq
for fastq in $IN_DIR/*.fq
do
  bn=`basename $fastq .fq`
  echo bwa_alignse_onlymapped.sh $WG_BWA_DB $fastq $OUT_DIR/$bn.bam
done
