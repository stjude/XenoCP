#!/bin/bash
# Performs STAR to align to a reference, converts to bam, and
# filters out all records for unmapped reads.
#
# $1 = reference db prefix
# $2 = input fastq file
# $3 = output bam

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOMEDIR=$1
FASTQ=$2
BAM=$3

# Write diagnostics
echo $0
hostname
date

# Make sure output directory exists
OUT_DIR=`dirname $BAM`
if [ "$OUT_DIR" != "" -a ! -e "$OUT_DIR" ]; then mkdir -p $OUT_DIR; fi

# Create temp dir to hold intermediate files
SCRATCH=`mktemp -d` || ( echo "Could not get scratch dir" >&2 ; exit 1 )

# Output filename
name=$(basename $BAM ".bam")

# Ensure reads are not compressed
reads=$(readlink -f ${FASTQ})
if [ $(file ${reads} | grep -c "gzip compressed data") -eq 1 ]
then
    gzip -dc ${reads} > reads.fq
    reads="reads.fq"
elif [ $(file ${reads} | grep -c "bzip2 compressed data") -eq 1 ]
then
    bzip2 -dc ${reads} > reads.fq
    reads="reads.fq"
fi

# STAR
n_cores=$(nproc)
cmd="STAR --genomeDir ${GENOMEDIR} \
             --readFilesIn ${reads} \
             --runMode alignReads \
             --runThreadN ${n_cores} \
             --outSAMunmapped Within \
             --outSAMstrandField intronMotif \
             --outSAMtype BAM Unsorted \
             --outSAMattributes NH HI AS nM NM MD XS \
             --outFilterMultimapScoreRange 1 \
             --outFilterMultimapNmax 20 \
             --outFilterMismatchNmax 10 \
             --alignIntronMax 500000 \
             --alignMatesGapMax 1000000 \
             --sjdbScore 2 \
             --alignSJDBoverhangMin 1 \
             --outFilterMatchNminOverLread 0.66 \
             --outFilterScoreMinOverLread 0.66 \
             --outFileNamePrefix ${name}. \
             --twopassMode Basic \
             --limitBAMsortRAM 3000000000 \
            && java.sh org.stjude.compbio.sam.TweakSam -V SILENT -G 4 -o ${BAM} -i ${name}.Aligned.out.bam"

echo $cmd 
if ! eval "$cmd"
then echo "STAR | TweakSam failed" >&2 ; exit 1
fi

# Remove scratch dir
rm -rf $SCRATCH

# Write end date
date
