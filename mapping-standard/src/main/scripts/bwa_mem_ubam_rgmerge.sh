#!/bin/bash
# Performs bwa mem to align to a reference, using unaligned BAM internally
# converted to FASTQ as input and merging the result with the unaligned BAM at
# the end; output is always coordinate sorted.
#
# The intermediate files are written to temp space.
#
# $1 = reference db prefix
# $2 = reference fasta file (must have a .dict equivalent file present as well)
# $3 = input ubam
# $4 = output bam
# $5 = "precopy" if you want the unaligned bam copied to scratch space first
# $6 = input FASTQ 1 if exists
# $7 = input FASTQ 2 if exists and input is paired

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PREFIX=$1
FASTA=$2
UBAM=$3
BAM=$4
PRECOPY=$5
FASTQ1=$6
FASTQ2=$7

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $BAM`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Create temp dir to hold intermediate files
SCRATCH=`mktemp -d` || ( echo "Could not get scratch dir" >&2 ; exit 1 )

# Pre-copy if specified
if [ "$PRECOPY" == "precopy" ]
then
  cp $UBAM $SCRATCH/ubam.bam
  UBAM=$SCRATCH/ubam.bam
fi

# Determine if we are doing pe or se mapping
paired=`samtools view $UBAM | awk 'and($2, 1) == 1 { print 1 } { exit }'`

# Generate FASTQs if needed
if [ "$FASTQ1" == "" ]
then
  FASTQ1=$SCRATCH/fq1.fq
  if [ $paired ]
  then
    FASTQ2="$SCRATCH/fq2.fq"
    fq2arg="F2=$FASTQ2"
  else
    FASTQ2=
    fq2arg=
  fi
  set -e -x
  java.sh picard.cmdline.PicardCommandLine SamToFastq INPUT=$UBAM FASTQ=$FASTQ1 $fq2arg RC=true
  set +e +x
# Check for paired ubam with unpaired FASTQ
elif [ "$FASTQ2" == "" -a $paired ]
then
  echo "You specified a single FASTQ file for a paired UBAM.  Please specify both FASTQs (or neither to auto-generate)" >&2
  exit 1
fi

# Intermediate aligned files
alignedsam=$SCRATCH/_aligned.sam
aligned=$SCRATCH/_aligned.bam

# Do the real work
set -e -x
bwa mem -M $PREFIX $FASTQ1 $FASTQ2 > $alignedsam
java.sh org.stjude.compbio.mapping.standard.CleanMappedSam -i $alignedsam -o $aligned -V SILENT -f $FASTA -t -e -d
java.sh org.stjude.compbio.mapping.standard.SimpleMergeBamAlignment -u $UBAM -a $aligned -o $BAM -d $FASTA.dict -P bwa_mem -O coordinate
set +e +x

# Write end date
date
