#!/bin/bash
# Tests a bam for existence, sorting/indexing, and stats available via flagstat
#
# Valid status values are:
# "unsorted" or "": do no status checking
# "queryname": sorted by queryname
# "coordinate": sorted by coordinate
# "indexed": sorted by coordinate, also check for index file
#
# Keep in mind that most of the real checking is done via the flagstat file, so
# if you omit the flagstat file or specify one that is out of date or for a
# different bam, then you may easily get bad results.
#
# $1 = bam file
# $2 = status (see above)
# $3 = (optional) flagstat.txt file; if omitted, no checks are done

# Source the test infrastructure script
source qclib.sh

# Process arguments
BAM=$1
STATUS=$2
FLAGSTAT=$3

# Start test case
starttestcase BamFile Bam: $BAM Status: $STATUS Flagstat: $FLAGSTAT

# Do tests

starttest BamExists
if [ ! -s "$BAM" ]; then aborttestcase; else passtest; fi

if [ -n "$FLAGSTAT" ]
then
  starttest FlagstatExists
  if [ ! -s "$FLAGSTAT" ]; then aborttestcase; else passtest; fi
fi

# Gets the declared sort order
# $1 = bam file
# Output: sort order: queryname, coordinate, unsorted, or ""
function getsorting {
  samtools view -H $1 | grep '@HD' | grep -o 'SO:[a-z]*' | sed 's/SO://'
}

if [ "$STATUS" == "queryname" ]
then
  starttest SortedByQueryname
  declsortorder=`getsorting $BAM`
  if [ "$declsortorder" == "queryname" ];
  then passtest;
  else failtest Got sort order [$declsortorder]
  fi
elif [ "$STATUS" == "coordinate" ]
then
  starttest SortedByCoordinate
  declsortorder=`getsorting $BAM`
  if [ "$declsortorder" == "coordinate" ];
  then passtest;
  else failtest Got sort order [$declsortorder]
  fi
elif [ "$STATUS" == "indexed" ]
then
  starttest SortedByCoordinateAndIndexed
  declsortorder=`getsorting $BAM`
  if [ "$declsortorder" == "coordinate" ];
  then
    if [ -f "$BAM.bai" ];
    then passtest; 
    else failtest No index file $BAM.bai;
    fi
  else failtest Got sort order [$declsortorder]
  fi
elif [ "$STATUS" != "" -a "$STATUS" != "unsorted" ]
then
  starttest ValidStatusSpecified
  failtest Invalid status specified: $STATUS
fi

if [ -s "$FLAGSTAT" ]
then
  starttest TotalReadCount
  readtot=`head -n 1 $FLAGSTAT | cut -d ' ' -f 1`
  if [ "$readtot" -gt 0 ]; then passtest; else failtest Got $readtot reads; fi
  
  starttest PairedReadCount
  paired=`grep 'paired in sequencing' $FLAGSTAT | cut -d ' ' -f 1`
  if [ "$paired" -eq 0 -o "$paired" -eq "$readtot" ]
  then passtest
  else failtest total=$readtot paired=$paired
  fi
  
  starttest Read1And2Counts
  read1=`grep -m 1 read1 $FLAGSTAT | cut -d ' ' -f 1`
  read2=`grep -m 1 read2 $FLAGSTAT | cut -d ' ' -f 1`
  if [ "$read1" -eq "$read2" ]
  then passtest
  else failtest read1=$read1 read2=$read2
  fi
  
  starttest Read1Plus2EqualsPaired
  let "read1and2 = read1 + read2"
  if [ "$paired" -eq "$read1and2" ]
  then passtest
  else failtest read1+read2=$read1and2 total=$readtot
  fi
  
  starttest MappingRate
  # Get the first whole number after an open paren
  mappingrate=`grep -m 1 mapped $FLAGSTAT | sed -r 's/^.*[(]([0-9]*).*$/\1/'`
  if [ "$mappingrate" -lt 65 ]
  then warntest Mapping rate $mappingrate% is less than 75%
  elif [ "$mappingrate" -lt 100 ]
  then passtest
  else warntest Mapping rate $mappingrate% is 100%
  fi
fi

summarize
