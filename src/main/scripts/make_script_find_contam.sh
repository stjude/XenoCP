#!/bin/bash
# Writes a script to create a list of reads from a bam file that are mapped
# better in at least one of a set of contamination bam files.
#
# Contamination BAMs must:
# * have the same basename as the matching input bam
# * have a subset of the reads of the matching input bam
# * be ordered in the same way as the matching input bam
#
# The output of the commands will be a list of read names called as
# contamination, and a bam of records that were not called contamination but
# which "tied" a contaminated record for manual review.  There is one list and
# one tie bam per input bam.
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = input directory containing sam/bams
# $2 = output directory for contamination lists
# $3 = output directory for tie bams
# $4... = input directories containing contamination sam/bams

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
IN_DIR=$1
OUT_DIR=$2
TIE_DIR=$3
shift 3
CONTAM_DIRS=$*

# Make sure output directory exists
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi
if [ "$TIE_DIR" != "" -a ! -e "$TIE_DIR" ]; then mkdir -p $TIE_DIR; fi

# One job per bam
for bam in $IN_DIR/*.bam
do
  bn=`basename $bam .bam`
  if [ "${bn:0:5}" == "noref" ]; then continue; fi
  contams=""
  for dir in $CONTAM_DIRS; do contams="$contams $dir/$bn.bam"; done
  echo java.sh org.stjude.compbio.xenocp.CreateContamLists -V SILENT -i $bam -o $OUT_DIR/$bn.contam.txt -t $TIE_DIR/$bn.tie.bam $contams
done
