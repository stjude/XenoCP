#!/usr/bin/env awk -f
# Converts SAM input data to FASTQ output, interleaving paired records in the
# same file.
# Pass nosuffix=1 to NOT add the /1 or /2 for PE reads (will mean duplicate
# read names)
/^[^@]/ {
  read = $1
  if(!nosuffix) {
    if(and($2, 0x40)) read = read "/1"
    else if(and($2, 0x80)) read = read "/2"
  }
  print "@" read
  print $10
  print "+"
  print $11  
}
