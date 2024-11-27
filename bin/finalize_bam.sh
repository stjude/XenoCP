#!/usr/bin/env bash
# Clean duplicate BAM header records from STAR aligned BAMs
# Add @PG record for XenoCP
#
# $1 = input bam
# $2 = output bam
# $3 = reference database
# $4 = aligner [star, bwa mem, bwa aln] used to search host genome
# $5 = original input BAM

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

BAM=$1
OUTPUT=$2

REFERENCE=$3
ALIGNER=$4
INPUT_BAM=$5

# Get BAM header and remove any consecutive duplicate records
# STAR aligned BAMs have a @CO record that gets duplicated 
# during XenoCP
samtools view -H ${BAM} | uniq > header.txt

# Get the ID value from the last @PG record
PG=$(awk '/^@PG/ {l=$0} END{match(l, /.*ID:([^\t]+).*/, arr); print arr[1] }' header.txt)

# Add a XenoCP @PG record to the BAM header
echo -e "@PG\tID:XenoCP\tPN:XenoCP\tCL:INPUT_BAM=${INPUT_BAM} \
ALIGNER=${ALIGNER} REFERENCE=${REFERENCE}\tDS:Unmap reads that \
are identified as deriving from the contaminant (host) genome \
\tPP:${PG}\tVN:4.0.0-alpha" >> header.txt

# Write the new header to an output BAM
samtools reheader header.txt ${BAM} > ${OUTPUT}

samtools index ${OUTPUT} ${OUTPUT}.bai

# Compute MD5 sum
if which md5sum > /dev/null
then
  md5sum ${OUTPUT} > ${OUTPUT}.md5
elif which md5 > /dev/null
then
  md5 -q ${OUTPUT} > ${OUTPUT}.md5
else
  echo "No MD5 program found" >&2; exit 1
fi