#!/bin/bash
# Writes a script to split an indexed, coordinate-sorted BAM with paired-end
# reads into pieces where both reads for a pair are guaranteed to be in the
# same file.  There will be multiple output files per chromosome.
#
# There will be a set of output files for each reference sequence in the bam 
# file, plus one prefixed noref and one prefixed crossref.  The named files will
# include all records for read pairs where both reads map to the same reference
# sequence.  The noref bams will contain unmapped reads with no reference, and
# the crossref bams will contain all reads whose mates are mapped to a different
# reference.
#
# The crossref bams take the longest to build.
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = input bam
# $2 = output directory
# $3 = number of pieces into which to split the crossref bam (default 5)
# $4 = number of pieces into which to split the noref bam (default 10)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
INPUT=$1
OUT_DIR=$2
CROSSREF_BUCKETS=$3
NOREF_BUCKETS=$4

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Apply defaults
if [ "$CROSSREF_BUCKETS" == "" ]; then CROSSREF_BUCKETS=5; fi
if [ "$NOREF_BUCKETS" == "" ]; then NOREF_BUCKETS=10; fi

# One job for cross-refs
echo tweak_split.sh $INPUT $OUT_DIR/crossref- bam -x -V SILENT -- -m -b $CROSSREF_BUCKETS -V SILENT

# One job for no-refs
echo tweak_split.sh $INPUT $OUT_DIR/noref- bam -X -n -V SILENT -- -m -b $NOREF_BUCKETS -V SILENT

# One job per sequence
seqs=`samtools view -H $INPUT | grep '^@SQ' | sed 's/^.*SN://' | cut -f 1`
for seq in $seqs
do echo tweak_split.sh $INPUT $OUT_DIR/ref- bam -X -c $seq -V SILENT -- -c -m -l 50000000 -b $NOREF_BUCKETS -V SILENT
done