#!/usr/bin/env bash
# A wapper of samtools view, sam_to_single.awk and picard SamToFastq for XenoCP pipeline
#
# Parameters:
# $1 = input aligned bam file
# $2 = output fastq file
#
# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 0; fi

inputBam=$1
outputFastq=$2

samtools view -h -F 4 $1 | sam_to_single.awk -v delim=. | java.sh picard.cmdline.PicardCommandLine SamToFastq INPUT=/dev/stdin FASTQ=$2 VALIDATION_STRINGENCY=SILENT QUIET=true VERBOSITY=ERROR
