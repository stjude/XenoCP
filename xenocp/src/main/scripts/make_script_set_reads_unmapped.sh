#!/bin/bash
# Writes a script to set reads with read names in a file to unmapped.
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = input bam
# $2 = output directory
# $3 = number of pieces into which to split the crossref bam (default 5)
# $4 = number of pieces into which to split the noref bam (default 10)

# Get parameters
INPUT_DIR=$1
NAMEFILES_DIR=$2
OUT_DIR=$3

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# One job per input dir bam
for bam in $INPUT_DIR/*.bam
do
  bn=`basename $bam .bam`
  namefile=`ls $NAMEFILES_DIR/$bn* 2> /dev/null`
  out=$OUT_DIR/$bn.bam
  if [ -f "$namefile" ]
  then
    echo java-settmp.sh - org.stjude.compbio.sam.TweakSam -i $bam -o $out -u $namefile -O coordinate -V SILENT
  else
    echo "No namefile for $bam.  Will write command to make symlink instead." >&2
    echo ln -s $bam $out
  fi
done

