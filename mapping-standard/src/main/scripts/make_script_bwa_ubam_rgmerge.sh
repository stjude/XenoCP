#!/bin/bash
# Writes a script to run bwa_alignxe_ubam_rgmerge.sh on a set of unaligned bam
# files in a directory; output is always coordinate sorted.
#
# Whether to use single- or paired-end mode is determined from the flags of the
# first record in each BAM
#
# Note: if the input filenames start with s_, then the s_ will be stripped when
# determining output filenames.
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = name of the configuration for the genome
# $2 = input directory containing unaligned bams
# $3 = output directory for bams
# $4 = "precopy" if you want the unaligned bams copied to tmp space first
# $5 = bwa alignment algorithm to use [backtrack or mem]
# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
CONFIG=$1
IN_DIR=$2
OUT_DIR=$3
ARGS=$@
PRECOPY=
ALGO="backtrack"
for arg in $ARGS
do
  if [ "$arg" == "precopy" ]; then PRECOPY=$arg
  elif [ "$arg" == "mem" ]; then ALGO="mem"
  elif [ "$arg" == "backtrack" ]; then ALGO="backtrack"
  fi
done

# Import configuration variables for whole genome bwa database and FASTA file
. import_config.sh genome $CONFIG WG_BWA_DB FASTA

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# One job per unaligned bam
for ubam in $IN_DIR/*.bam
do
  # Determine output location and check for existing file
  bam=$OUT_DIR/`basename $ubam | sed 's/^s_//'`
  if [ -f "$bam" ]; then echo "Skipping existing bam: $bam" >&2; continue; fi
  
  # Write the command
  if [ "$ALGO" == "backtrack" ] 
  then 
    echo bwa_alignxe_ubam_rgmerge.sh $WG_BWA_DB $FASTA $ubam $bam $PRECOPY
  elif [ "$ALGO" == "mem" ]
    then
      echo bwa_mem_ubam_rgmerge.sh $WG_BWA_DB $FASTA $ubam $bam $PRECOPY
  fi
done
