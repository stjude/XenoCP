#!/bin/bash
# Writes a script to run merge_markdup_index.sh on a set of aligned bam
# files in a directory using a run config file as a guide
#
# The output files will be named SAMPLE.bam or SAMPLE-BUILD.bam.  BUILD will
# be included if the build in the run config is not empty, "-", ".", or
# "LATEST".
#
# The commands are written to stdout, so they must be redirected to a file if
# you want to write a script file.
#
# $1 = run configuration file
# $2 = root directory containing sample directories
# $3 = subpath to aligned BAMs dir, relative to sample dir
# $4 = subpath to output dir, relative to sample dir

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
RUN_CONFIG=$1
ROOT_DIR=$2
IN_SUBPATH=$3
OUT_SUBPATH=$4

# Loop through run config entries
cat $RUN_CONFIG | while read sample indexes genome build bam ignore
do
  # Transform build
  build=`echo $build | sed 's/^[-.]$//'`
  
  # Get directories
  indir=$ROOT_DIR/$sample/$IN_SUBPATH
  outdir=$ROOT_DIR/$sample/$OUT_SUBPATH
  
  # Creat output directory if it doesn't exist
  mkdir -p $outdir
  
  # Make index list into a character set for globbing
  indexset=`echo $indexes | tr -d ,`
  
  # Make output BAM filename
  outname=$sample
  if [ "$build" != "" -a "$build" != "LATEST" ]
  then outname=$outname-$build
  fi
  outname=$outname.bam

  if [ "$bam" == "-" ]
  then 
   bam=
  fi
 
  # Write the MMI command
  echo merge_markdup_index.sh $outdir/$outname $indir/[$indexset]*.bam $bam
done
