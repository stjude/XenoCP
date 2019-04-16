#!/bin/awk -f
# Removes pairing information from SAM data.
#
# Read 1/2 information is added as a suffix on the read name.
#
# Parameters (all are optional):
# - nosuffix: pass true value to NOT add the /1 or /2 for PE reads (will mean 
#   duplicate read names)
# - delim: pass an alternate delimiter to use in place of /
BEGIN {
  FS = "\t"
  OFS = FS
  if(!delim) delim = "/"
}
/^@/ {
  print
  next
}
{
  read = $1
  if(!nosuffix) {
    if(and($2, 0x40)) read = read delim "1"
    else if(and($2, 0x80)) read = read delim "2"
  }
  # pairing-related: 1, 2, 8, 20, 40, 80
  flags = and($2, compl(0xEB))
  
  # Reset and print
  $1 = read
  $2 = flags
  $7 = "*"
  $8 = "0"
  $9 = "0"
  print
}