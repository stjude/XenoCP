#!/bin/sh
# Compiles mapping stats by reading flagstat files, and writes them to stdout.
#
# The arguments are either flagstat.txt files, or sample/directory names.  If
# non-files are given, then the flagstat.txt file is searched for in the
# following locations (in order), based on the specified argument ARG:
# 1. ARG.flagstat.txt
# 2. ARG/flagstat.txt
# 3. ARG/xcp/bam-purified/*flagstat.txt
# 4. ARG/finish/*flagstat.txt
# 5. ARG/resolve-merged/*flagstat.txt
#
# The first location where a file is found is chosen
#
# $1 = (optional) -h for fixed-width, human-readable output
# $1,... = flagstat.txt files or directories containing them

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
HR=$1
if [ "$HR" == "-h" ]; then shift; else HR= ; fi
INPUTS="$@"

if [ $HR ]
then printf "%-40s %13s %13s %13s  %5s  %5s\n" Sample Reads Mapped NonDupMapped Mpd% Dup%
else echo -e "Sample\tReads\tMapped\tNonDupMpd\tMpdRatio\tDupRatio"
fi
for input in $INPUTS
do
  # Find flagstat file(s)
  if [ -f "$input" ]; then fslist=$input
  elif [ -f "$input.flagstat.txt" ]; then fslist=$input.flagstat.txt
  elif [ -f "$input/flagstat.txt" ]; then fslist=$input/flagstat.txt
  elif [ "`ls $input/xcp/bam-purified/*flagstat.txt 2> /dev/null`" ]; then fslist=`ls $input/xcp/bam-purified/*flagstat.txt`
  elif [ "`ls $input/finish/*flagstat.txt 2> /dev/null`" ]; then fslist=`ls $input/finish/*flagstat.txt`
  elif [ "`ls $input/resolve-merged/flagstat.txt 2>/dev/null`" ]; then fslist=`ls $input/resolve-merged/*flagstat.txt`
  else echo "Skipping $input" >&2; continue
  fi
  
  # Loop over them and print info
  for fs in $fslist
  do
    # Determine the label (filename w/out .flagstat.txt suffix if possible,
    # otherwise, the basename of the dirname)
    label=`basename $fs .flagstat.txt`
    if [ "$label" == "flagstat.txt" ]
    then
      label=`dirname $fs`
      label=`basename $label`
    fi

    # Works with both versions of flagstat output
    awk -v h=$HR -v label=$label '
      BEGIN { OFS="\t" }
      /total/ { total=$1 } /duplicates/ {  dup=$1; getline; mapped=$1 } 
      END {
        nondupmapped = mapped - dup
        if(total == 0) pct = 0; else pct = mapped / total
        if(mapped == 0) duppct = 0; else duppct = dup / mapped
        if(h == "") {
          print label, total, mapped, nondupmapped, pct, duppct
        }
        else {
          pct = 100 * pct
          duppct = 100 * duppct
          printf("%-40s %13'"'"'d %13'"'"'d %13'"'"'d  %4.1f%%  %4.1f%%\n", label, total, mapped, nondupmapped, pct, duppct)
        }
      }' $fs
  done 
done
