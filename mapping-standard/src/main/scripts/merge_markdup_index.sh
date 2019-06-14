#!/usr/bin/env bash
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
# $1 = (optional) -t <threads> to execute steps multithreaded
# $2 = output bam/sam
# $3... = input bam/sam files

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

NUM_THREADS=1
NO_MARKDUP=

# Get parameters
while (( "$#" )); do
  case "$1" in
    --no-markdup)
      NO_MARKDUP=$1
      shift
      ;;
    -t|--threads)
      NUM_THREADS=$2
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

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
cmd="sambamba merge -t $NUM_THREADS $L_MERGED $INPUTS"
echo $cmd
if ! $cmd
then echo "Merge command failed" >&2; exit 1
fi

# Do the mark duplicates/index/md5 (cannot pipe this because MarkDuplicates makes multiple input passes)
if [ $NO_MARKDUP ]
then
  echo "Skipping MarkDuplicates because $NO_MARKDUP was specified"
else
  cmd="sambamba markdup -t $NUM_THREADS $L_MERGED $L_OUTPUT"
  echo $cmd
  if ! $cmd
  then echo "Mark duplicates command failed" >&2; exit 1
  fi
fi

sambamba index -t $NUM_THREADS $L_OUTPUT ${bn}.bai
sambamba flagstat -t $NUM_THREADS $L_OUTPUT  > $FLAGSTAT
if which md5sum > /dev/null
then
  md5sum $L_OUTPUT > $L_OUTPUT.md5
elif which md5 > /dev/null
then
  md5 -q $L_OUTPUT > $L_OUTPUT.md5
else
  echo "No MD5 program found" >&2; exit 1
fi

# Copy output
cmd="cp $SCRATCH_DIR/$bn.bam $SCRATCH_DIR/$bn.bam.bai $SCRATCH_DIR/$bn.bam.md5 $OUT_DIR/"
echo $cmd
if ! $cmd
then echo "Copy to output failed" >&2; exit 1
fi
