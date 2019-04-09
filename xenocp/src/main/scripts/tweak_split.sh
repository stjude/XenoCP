#!/bin/bash
# Calls TweakSam and pipes into SplitSam; useful for parallized BAM splitting
#
# Example call to extract read pairs both mapped to chr 1 and then smartly split
# them into 10 buckets:
#
# tweak_split.sh input.bam outdir/output bam -X -c 1 -- -m -b 10
#
# $1 = input sam/bam file
# $2 = output sam/bam prefix with path and trailing delimiter (not a bare dir)
# $3 = output sam/bam extension, without dot (this determines format)
# $4,... = TweakSam parameters
# -- = use literal -- to delineate TweakSam parameters from SplitSam parameters
# $N,... = SplitSam parameters

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
INPUT=$1
OUT_PREFIX=$2
OUT_EXT=$3
shift 3
TWEAK_PARAMS=
SPLIT_PARAMS=
ONSPLIT=
for PARAM in $*
do
  if [ "$PARAM" == "--" ]
  then
    ONSPLIT=1
    continue
  fi
  if [ -z "$ONSPLIT" ]
  then
    TWEAK_PARAMS="$TWEAK_PARAMS $PARAM"
  else
    SPLIT_PARAMS="$SPLIT_PARAMS $PARAM"
  fi
done

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $OUT_PREFIX`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Do the work
tcmd="java.sh org.stjude.compbio.sam.TweakSam -i $INPUT $TWEAK_PARAMS"
scmd="java.sh org.stjude.compbio.sam.SplitSam $SPLIT_PARAMS $OUT_PREFIX $OUT_EXT"
echo $tcmd \| $scmd
$tcmd | $scmd
