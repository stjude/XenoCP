#!/bin/bash

BAM=$1

samtools view -H ${BAM} | grep @SQ | cut -f 2  | cut -f 2 -d':' 
#for chr in $(samtools view -H ${BAM} | grep @SQ | cut -f 2  | cut -f 2 -d':')
#do
#echo -n "$chr,"
#done

#echo -n  "{\"chroms\": [";

#for chr in $(samtools view -H ${BAM} | grep @SQ | cut -f 2  | cut -f 2 -d':')
#do
#  echo -n "\"$chr\","
#done

#echo -n "]}"

