#!/bin/sh
# Performs merge, mark duplicates, index of a set of sorted bam files using
# Picard.  Also calculates flagstat and md5.
#
# IMPORTANT: All inputs must be coordinate sorted, although the headers do
# not need to declare it.
#
# Input is NOT copied to scratch.  Scratch is used for Picard temp files and
# for output.
#
# $0 = (optional) --no-markdup to skip the mark dup step
# $1 = output bam/sam
# $2... = input bam/sam files

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
if [ "$1" == "--no-markdup" ]
then NO_MARKDUP=$1; shift
else NO_MARKDUP=
fi
OUTPUT=$1
shift
INPUTS=$*

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $OUTPUT`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Get conventionally-named output locations
bn=`basename $OUTPUT .bam`
FLAGSTAT=$OUT_DIR/$bn.flagstat.txt
MARKDUP_METRICS=$OUT_DIR/$bn.markdup.txt

# Create scratch dir
SCRATCH_DIR=`mktemp -d`

# Get local output locations
L_OUTPUT=$SCRATCH_DIR/$bn.bam
if [ $NO_MARKDUP ]
then L_MERGED=$L_OUTPUT
else L_MERGED=$SCRATCH_DIR/merged.bam
fi

# Do the merge
INPUT_ARGS=
for INPUT in $INPUTS; do INPUT_ARGS="$INPUT_ARGS INPUT=$INPUT"; done
cmd="java-settmp.sh $SCRATCH_DIR picard.cmdline.PicardCommandLine MergeSamFiles VALIDATION_STRINGENCY=LENIENT $INPUT_ARGS OUTPUT=$L_MERGED ASSUME_SORTED=true CREATE_INDEX=true CREATE_MD5_FILE=true"
date
echo $cmd
if ! $cmd
then echo "Merge command failed" >&2; exit 1
fi

# Do the mark duplicates/index/md5 (cannot pipe this because MarkDuplicates makes multiple input passes)
if [ $NO_MARKDUP ]
then
  echo "Skipping MarkDuplicates because $NO_MARKDUP was specified"
else
  cmd="java-settmp.sh $SCRATCH_DIR -Xmx8g picard.cmdline.PicardCommandLine MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=$L_MERGED OUTPUT=$L_OUTPUT METRICS_FILE=$MARKDUP_METRICS ASSUME_SORTED=true CREATE_INDEX=true CREATE_MD5_FILE=true"
  date
  echo $cmd
  if ! $cmd
  then echo "Mark duplicates command failed" >&2; exit 1
  fi
fi

# Copy output
cmd="cp $SCRATCH_DIR/$bn.bam $SCRATCH_DIR/$bn.bai $SCRATCH_DIR/$bn.bam.md5 $OUT_DIR/"
date
echo $cmd
if ! $cmd
then echo "Copy to output failed" >&2; exit 1
fi

# Rename the index
mv $OUT_DIR/$bn.bai $OUT_DIR/$bn.bam.bai

# Do the flagstat
cmd="samtools flagstat $L_OUTPUT"
date
echo $cmd \> $FLAGSTAT
if ! $cmd > $FLAGSTAT
then echo "Warning: flagstat command failed" >&2
fi
