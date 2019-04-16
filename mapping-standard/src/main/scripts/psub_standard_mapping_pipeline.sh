#!/bin/bash
# Submits a standard mapping job to Pallas. 
#
# $0 = (optional) --no-run to only create and not run the pipeline in pallas
# $0 = (optional) --no-psub to only write the psub command and not run it
# $1 = Absolute path to run directory
# $2 = Pallas queue (e.g. lsf-pcgp-low, not pcgp_low)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
while [ "${1:0:2}" == "--" ]
do
  if [ "$1" == "--no-run" ]
  then noRun=$1
  elif [ "$1" == "--no-psub" ]
  then noPsub=$1
  else echo "Unrecognized option: $1" >&2; exit 1
  fi
  shift
done
run_dir=$1
queue=$2

# Get run name from run directory
run_name=`basename $run_dir`
# Use the absolute path to get the relative path (used by pallas)
run_dir_relative_path=`echo $run_dir | sed 's#/nfs_exports/genomes/1/##'`
# User name
user=`whoami`
# Build full json string (201457594 is the id of the std anls template)
json=`tr -d '\n' << EOF
{
"p-template-id":"201457594",
"p-user-name":"$user",
"p-project-name":"$LSB_DEFAULTPROJECT",
"p-requested-queue-name":"$queue",
"p-pipeline-name":"Std Mapping for $run_name",
"run_name":"$run_name",
"vfs_dir":"genomes/$run_dir_relative_path",
"run_dir":"$run_dir",
"PATH":"$PATH",
"CLASSPATH":"$CLASSPATH",
"PERL5LIB":"$PERL5LIB",
"R_LIBS":"$R_LIBS",
"SJ_CONFIGS":"$SJ_CONFIGS",
"script01c":"$run_dir/01c.sh",
"cmds01":"$run_dir/cmds-01.sh",
"script02a":"$run_dir/02a.sh",
"script02c":"$run_dir/02c.sh",
"cmds02":"$run_dir/cmds-02.sh",
"script03a":"$run_dir/03a.sh",
"script03b":"$run_dir/03b.sh"
}
EOF`

# Build the URL
BASE_URL=http://pallas.stjude.org/service/pipeline
if [ $noRun ]
then action=createPipeline
else action=createAndRunPipeline
fi
url=$BASE_URL/$action

# Perform the curl call
cmd="curl --include --request POST --header 'Content-Type:application/json' $url -d '$json'"
echo Command: "$cmd"
if [ $noPsub ]
then
  echo "psub was disabled via $noPsub, nothing is submitted to pallas"
else 
  eval "$cmd"
fi
