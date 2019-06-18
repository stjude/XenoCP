#!/bin/bash
# Performs Picard SortSam and then flagstat
#
# Flagstat file will have same name as output BAM, but with .flagstat.txt in
# place of .bam
#
# Parameters:
# $1 = input file
# $2 = sort order
# $3 = output BAM file
# $4,... = other parameters are sent to SortSam

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
INPUT=$1
SORT_ORDER=$2
OUTPUT=$3
shift 3
ARGS="$@"

# Make output dir if it does not exist
OUT_DIR=`dirname $OUTPUT`
if [ ! -d "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Run SortSam
cmd="java.sh picard.cmdline.PicardCommandLine SortSam INPUT=$INPUT OUTPUT=$OUTPUT SORT_ORDER=$SORT_ORDER VALIDATION_STRINGENCY=SILENT $ARGS"
echo "$cmd"
if ! eval "$cmd"
then echo "Error while running SortSam" >&2; exit 1
fi

# Run flagstat
flagstat_file=$OUT_DIR/`basename $OUTPUT .bam`.flagstat.txt
cmd="samtools flagstat $OUTPUT > $flagstat_file"
echo "$cmd"
eval "$cmd"