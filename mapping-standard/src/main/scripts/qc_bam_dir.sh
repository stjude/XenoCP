#!/bin/bash
# Tests a dir that should contain one output BAM for every input BAM or pair of
# FASTQs
#
# $1 = directory that contains output bam files
# $2 = directory that contains input ubam or fastq files
# $3 = (optional) required sorting: coordinate or queryname

# Source the test infrastructure script
source qclib.sh

# Process arguments
DIR=$1
SEQDIR=`readlink -f $2`
SORTORDER=$3

# Start test case
starttestcase OneTilegrpBamsDir Dir: $DIR Seq Dir: $SEQDIR

# Do tests

starttest DirExists
if [ ! -d "$DIR" ]; then aborttestcase; else passtest; fi

starttest SeqDirExists
if [ ! -d "$SEQDIR" ]; then aborttestcase; else passtest; fi

cd $DIR

starttest SomeBamFilesExist
if ls *.bam >&2
then passtest
else aborttestcase No output bam files
fi

starttest SeqDirContainsFiles
cnt=`ls $SEQDIR/s_?_1_*.fq 2>/dev/null | wc -l`
SEQTYPE=
if [ "$cnt" -gt 0 ]
then
  SEQTYPE=PEFASTQ
  passtest
else
  cnt=`ls $SEQDIR/s_?_*.fq 2>/dev/null | wc -l`
  if [ "$cnt" -gt 0 ]
  then
    SEQTYPE=SEFASTQ
    passtest
  else
    cnt=`ls $SEQDIR/s_?_*.bam 2>/dev/null | wc -l`
    if [ "$cnt" -gt 0 ]
    then
      SEQTYPE=UBAM
      passtest
    else
      aborttestcase No fastq or ubam sequence files
    fi
  fi
fi

starttest AllBamFilesExist
bamcnt=`ls set_*.bam 2>/dev/null | wc -l`
BAMTYPE=s
if [ "$bamcnt" -eq "$cnt" ]
then passtest
elif [ "$bamcnt" -gt 0 ]
then failtest Expected $cnt bams but found $bamcnt
else
  BAMTYPE=n
  bamcnt=`ls *.bam 2>/dev/null | grep -cE '^._'`
  if [ "$bamcnt" -eq "$cnt" ]
  then passtest
  else failtest Expected $cnt bams but found $bamcnt
  fi
fi

# Will do bam file size tests one set at a time, so make a function to help
# with that:

# Tests the bam file sizes for one set.
# Uses find_outlier_files.pl to find file with a size that is too small, and
# calls them as bad, unless the corresponding seq file is also small
#
# Params:
# $1 = set number
function test_size_for_one_set {
  index=$1
  
  starttest BamFileSize-set_$index
  
  # Get bam and seq arrays
  if [ "$BAMTYPE" == "s" ]
  then bams=( set_${index}_*.bam )
  else bams=( ${index}_*.bam )
  fi
  if [ "$SEQTYPE" == "PEFASTQ" ]
  then seqs=( $SEQDIR/s_${index}_1_*.fq )
  elif [ "$SEQTYPE" == "SEFASTQ" ]
  then seqs=( $SEQDIR/s_${index}_1_*.fq )
  else seqs=( $SEQDIR/s_${index}_*.bam )
  fi
  count=${#bams[@]}
  
  # If there's only one BAM, then warn
  if [ "$count" == "1" ]
  then
    #warntest Only one BAM, probably too little yield
    echo Only one BAM, possibly too little yield >&2
    passtest
    return
  fi
  
  # Compute bam size/seq size (b/s) for each matching bam/seq, as well as the max
  bsx100=()
  maxbsx100=0
  i=0
  let "last = count - 1"
  while [ "$i" -lt "$count" ]; do
    # File sizes and ratio
    b=`filesize ${bams[$i]}`
    s=`filesize ${seqs[$i]}`
    bsx100[$i]=`expr 100 \* $b / $s`
    # Update max if a new max is found (except on last bam)
    if [ "$i" -lt "$last" -a "${bsx100[$i]}" -gt "$maxbsx100" ]
    then maxbsx100=${bsx100[$i]}
    fi
    # Increment i for loop
    let "++i"
  done

  # Gives an idea of how reliable this will be  
  echo "bsx100=${bsx100[@]} maxbsx100=$maxbsx100" >&2
  
  # Check for maxbsx100 = 0 (indicates no bams or all are very small)
  if [ "$maxbsx100" -lt 1 ]
  then
    failtest No BAMs or else all BAMs were very small
    return
  fi
  
  # Now, for each bam, use max b/s - b/s to make a call
  warn=
  i=0
  while [ "$i" -lt "$last" ]; do
    curbsx100=${bsx100[$i]}
    #let "dbsx100 = maxbsx100 - curbsx100"
    let "cbsx100 = 100 * curbsx100 / maxbsx100"
    #if [ "$dbsx100" -ge 25 ]
    if [ "$cbsx100" -lt 30 ]
    then
      failtestifactive At least one BAM appears to be too small
      echo "bam=${bams[$i]} fq=${fqs[$i]} 100*(bam size)/(fq size)=bsx100=$curbsx100 max(bsx100)=$maxbsx100" >&2
    #elif [ "$dbsx100" -ge 10 ]
    elif [ "$cbsx100" -lt 50 ]
    then
      warn=1
      echo "bam=${bams[$i]} fq=${fqs[$i]} 100*(bam size)/(fq size)=bsx100=$curbsx100 max(bsx100)=$maxbsx100" >&2
    fi
    let "++i"
  done
  if [ "${bsx100[$last]}" -lt 5 ]
  then
    failtestifactive Last BAM appears to be too small
  elif [ "${bsx100[$last]}" -lt 30 ]
  then
    warn=1
  fi
  if [ -z "$warn" ]
  then passtestbydefault
  else warntestbydefault At least one BAM size was questionable
  fi
}

# Do BAM size tests one set at a time
if [ "$BAMTYPE" == "s" ]
then indexes=`ls set_*.bam | cut -c 5 | sort | uniq`
else indexes=`ls [0-9,A-Z]*.bam | cut -c 1 | sort | uniq`
fi
for index in $indexes; do test_size_for_one_set $index; done

# List of all bams for looping
if [ "$BAMTYPE" == "s" ]
then bams=`ls set_*.bam`
else bams=`ls *.bam | grep -E '^[0-9,A-Z]'`
fi

if [ -n "$SORTORDER" ]
then
  starttest DeclaredSortOrder
  
  for bam in $bams
  do
    declsortorder=`samtools view -H $bam | grep '@HD' | grep -o 'SO:[a-z]*' | sed 's/SO://'`
    if [ "$declsortorder" != "$SORTORDER" ]
    then failtestifactive Declared sort order was [$declsortorder], not [$SORTORDER] as expected in $bam
    fi
  done
  
  passtestbydefault
fi

starttest MissingRefName
# Actually, will only check one (there's a break at the end)
for bam in $bams
do
  numbad=`samtools view -F 4 $bam | head | cut -f 3 | grep -c '*'`
  if [ "$numbad" -ne "0" ]; then failtest Found mapped records with no ref name
  else passtest
  fi
  break
done

summarize
