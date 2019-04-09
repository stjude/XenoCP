# Tests that a BAM was successfully split into component BAMs
#
# Note that this only tests read counts to try to ensure that all of the reads
# are accounted for.  It does not check that any other requirement of the
# splitting is fulfilled
#
# $1 = Directory full of split BAM's flagstats (each file must match *flagstat*)
# $2 = Input BAM's flagstat

# Source the test infrastructure script
source qclib.sh

# Process arguments
DIR=$1
INPUT=$2

# Start test case
starttestcase BamSplitResults Dir: $DIR Input: $INPUT

# Do tests

starttest DirExists
if [ ! -d "$DIR" ]; then aborttestcase; else passtest; fi
starttest InputFsExists
if [ ! -f "$INPUT" ]; then aborttestcase; else passtest; fi

# Do calculations into temp files
# Pipe flagstats in, writes summary out
fssum() {
  awk '{ word = $2; if(word == "+") word = $4; sum[word] += $1 } END { for(word in sum) print sum[word], word }'
}
insum=`mktemp`
cat $INPUT | fssum > $insum
outsum=`mktemp`
cat $DIR/*flagstat* | fssum > $outsum

# Check each of the counts that we care about
# Writes count1\ncount2
# $1 = category
getcounts() {
  echo `cat $insum $outsum | awk -v cat=$1 '$2 == cat { print $1 }'`
}

starttest TotalReads
read incount outcount < <(getcounts in)
if   [ "$incount" == "" ]; then failtest Count not get read count
elif [ "$incount" == "0" ]; then warntest No input reads
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Total read counts differ: input=$incount output=$outcount
fi

starttest MappedReads
read incount outcount < <(getcounts mapped)
if   [ "$incount" == "" ]; then failtest Count not get mapped read count
elif [ "$incount" == "0" ]; then warntest No input mapped reads
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Mapped read counts differ: input=$incount output=$outcount
fi

starttest DuplicateReads
read incount outcount < <(getcounts duplicates)
if   [ "$incount" == "" ]; then failtest Count not get duplicate read count
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Duplicate read counts differ: input=$incount output=$outcount
fi

starttest PairedReads
read incount outcount < <(getcounts paired)
if   [ "$incount" == "" ]; then failtest Count not get paired read count
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Paired read counts differ: input=$incount output=$outcount
fi

starttest Reads1
read incount outcount < <(getcounts read1)
if   [ "$incount" == "" ]; then failtest Count not get read1 read count
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Read 1 read counts differ: input=$incount output=$outcount
fi

starttest Reads2
read incount outcount < <(getcounts read2)
if   [ "$incount" == "" ]; then failtest Count not get read2 read count
elif [ "$incount" == "$outcount" ]; then passtest
else failtest Read 2 read counts differ: input=$incount output=$outcount
fi

rm $insum $outsum

summarize
