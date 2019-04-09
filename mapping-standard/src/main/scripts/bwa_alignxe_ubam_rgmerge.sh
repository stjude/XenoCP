#!/bin/bash
# Performs bwa aln and bwa sampe/-se to align to a reference, using unaligned
# BAM as input and merging the result with the unaligned BAM at the end; output
# is always coordinate sorted.
#
# The intermediate file is written to temp space.
#
# $1 = reference db prefix
# $2 = reference fasta file (must have a .dict equivalent file present as well)
# $3 = input ubam
# $4 = output bam
# $5 = "precopy" if you want the unaligned bam copied to scratch space first

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PREFIX=$1
FASTA=$2
UBAM=$3
BAM=$4
PRECOPY=$5

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $BAM`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Create temp dir to hold intermediate files
SCRATCH=`mktemp -d --tmpdir` || ( echo "Could not get scratch dir" >&2 ; exit 1 )

# Pre-copy if specified
if [ "$PRECOPY" == "precopy" ]
then
  cp $UBAM $SCRATCH/
  UBAM=$SCRATCH/`basename $UBAM`
fi

# Intermediate aligned files
alignedsam=$SCRATCH/_aligned.sam
aligned=$SCRATCH/_aligned.bam

# Determine if we are doing pe or se mapping
paired=`samtools view $UBAM | awk 'and($2, 1) == 1 { print 1 } { exit }'`

# Perform alignment
if [ $paired ]
then
  # bwa aln 1/2
  for M in 1 2
  do
    cmd="bwa aln -b -$M $PREFIX $UBAM"
    sai=$SCRATCH/`basename $BAM .bam`_$M.sai
    eval sai$M=$sai
    date
    echo $cmd \> $sai
    if ! $cmd > $sai
    then echo "Failure in bwa aln $M" >&2 ; exit 1
    fi
  done
  
  # bwa sampe
  bwacmd="bwa sampe $PREFIX $sai1 $sai2 $UBAM $UBAM"
  date
  echo "$bwacmd > $alignedsam"
  if ! $bwacmd > $alignedsam
  then echo "Failure in bwa sampe" >&2 ; exit 1
  fi
else
  # bwa aln
  cmd="bwa aln -b $PREFIX $UBAM"
  sai=$SCRATCH/`basename $BAM .bam`.sai
  date
  echo $cmd \> $sai
  if ! $cmd > $sai
  then echo "Failure in bwa aln" >&2 ; exit 1
  fi

  # bwa samse
  bwacmd="bwa samse $PREFIX $sai $UBAM"
  date
  echo "$bwacmd > $alignedsam"
  if ! $bwacmd > $alignedsam
  then echo "Failure in bwa samse" >&2 ; exit 1
  fi
fi

# CleanSam (and convert to BAM)
cmd="java.sh org.stjude.compbio.mapping.standard.CleanMappedSam -i $alignedsam -o $aligned -V LENIENT -f $FASTA -t -e -d -P"
date
echo $cmd
if ! $cmd
then echo "Failure in CleanMappedSam" >&2 ; exit 1
fi

# MergeBamAlignment
mbacmd="java.sh org.stjude.compbio.mapping.standard.SimpleMergeBamAlignment -u $UBAM -a $aligned -o $BAM -d $FASTA.dict -P bwa -O coordinate"
date
echo $mbacmd
if ! $mbacmd
then echo "Failure in SimpleMergeBamAlignment" >&2 ; exit 1
fi

# Remove scratch dir
rm -rf $SCRATCH
#DEBUG: echo $SCRATCH

# Write end date
date
