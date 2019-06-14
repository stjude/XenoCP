#!/bin/bash

BAM=$1

samtools view -H ${BAM} | grep @SQ | cut -f 2  | cut -f 2 -d':' 
