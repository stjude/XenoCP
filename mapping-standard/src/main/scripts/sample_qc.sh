#!/bin/bash
# The mapping or coverage statistics are piped into this script, along with a
# cutoff value.
#
# The input to be piped in is a 2-column whitespace-delimited table where the
# first column is the sample name and the second is the statistic.
#
# Each sample is printed to STDOUT with a PASS/FAIL indicator. 
#
# You must also specify the following parameters:
#
# $1 = name of the test (use quotes if nec.), e.g. "Mapping Stats" or Coverage
# $2 = name of the statistic, e.g. NonDupMapped or PctExonsGt20x
# $3 = cutoff value, e.g. 50000000 or 80
#
# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
TEST=$1
NAME=$2
CUTOFF=$3

# Get failed qc file
awk -v cutoff=$CUTOFF -v name=$NAME '
  BEGIN          { OFS="\t" }
  $1 == "Sample" { next }
  $2 < cutoff    { print "FAIL", $1, $2, cutoff }
  $2 >= cutoff   { print "PASS", $1, $2, cutoff }
  ' > sample_qc.txt

