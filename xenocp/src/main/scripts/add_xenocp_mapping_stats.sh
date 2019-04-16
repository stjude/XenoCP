#!/bin/sh
# Takes mapping stats from compile_mapping_stats.sh and appends the xenograft
# cleansing statistics.
#
# The mapping stats file must have been built using sample/directory names.  If
# it was built with full filenames, then this will fail.
#
# The sample names are read from the given mapping stats file, and a second set
# of flagstat files are searched for in these locations, in this order (where
# SAMPLE is the sample name read from the first column of the input file):
# 1. SAMPLE/refine/bam/*flagstat.txt
# 2. SAMPLE/WG/*flagstat.txt
# 3. ../BucketRaw/QC/SAMPLE.flagstat.txt
#
# The first location where a file is found is chosen.
#
# Output is written to stdout; be sure you do not redirect it over the original
# file, as that will not work (the original file will be overwritten before this
# script runs, and it will have no input)
#
# $1 = mapping stats file from compile_mapping_stats.sh

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
INPUT=$1

# Start temp file containing new columns
newcolfile=`mktemp`
echo -e "Removed (xeno)\tRemoved Ratio" > $newcolfile

# Build it up
cut -f 1,3 $INPUT | while read input mapped
do
  # Find the pre-cleansing flagstat
  if [ -f "`ls $input/refine/bam/*flagstat.txt 2> /dev/null`" ]; then fs=`ls $input/refine/bam/*flagstat.txt`
  elif [ -f "`ls $input/WG/*flagstat.txt 2>/dev/null`" ]; then fs=`ls $input/WG/*flagstat.txt`
  elif [ -f "`ls ../BucketRaw/QC/$input.flagstat.txt 2>/dev/null`" ]; then fs=`ls ../BucketRaw/QC/$input.flagstat.txt`
  else echo "Skipping $input" >&2; echo; continue
  fi

  # Write the original number mapped and the removal ratio to the temp file
  awk -v mapped=$mapped 'BEGIN { OFS="\t" } NR == 4 { removed = $1 - mapped; print removed, removed / mapped; exit 0 }' $fs >> $newcolfile
done

# Add the columns to the original file and write to stdout
paste $INPUT $newcolfile

# Remove temp file
rm $newcolfile
