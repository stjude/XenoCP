#!/bin/bash
# Performs several setup tasks for a run of the Standard Mapping Pipeline on
# LSF using a config file.
#
# The config file is tab-delimited text with the following columns:
# 1. Sample: sample name
# 2. Index list: comma-separated list of indexes in the range 0-9.
#    This tells which UBAMs to include in the final BAM
# 3. Genome config name (e.g. GRCh37-lite)
# 4. Build: build name to incorporate into the BAM filename; if you specify
#    an empty string, "-", ".", or "LATEST", then none is added
# Any additional fields are ignored, so you can put more fields for other use,
# if you would like.
#
# Note that all mapping will be done using the same genome config, which is 
# specified as a parameter.  So, column 3 is ignored in all rows.
#
# Also, the config file may contain comment lines starting with a #.  These will
# be ignored.
#
# $1 = run config file
# $2 = genome name
# $3 = threshold for mapping cutoff
# $4 = run directory
# $5 = data directory
# $6 = [optional] bwa algorithm to use

# Options
# -i input subpath
# -int intermediate subpath
# -o output subpath
# -r run subpath

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
RUN_CONFIG=$1
GENOME_NAME=$2
CUTOFF=$3
RUN_DIR=$4
DATA_DIR=$5
ALGO=$6
if [ "$ALGO" == "" ]; then ALGO="backtrack"; fi

# Get config name and run dir name
RUN_CONFIG_PATH=`readlink -f $RUN_CONFIG`
RUN_CONFIG_NAME=`basename $RUN_CONFIG`

if [ "$CUTOFF" == "" ] 
then 
  CUTOFF=0
fi

# Read pipeline config
. import_config.sh app mapping-standard

# Read step script
. steplib.sh

if [ ! -d $RUN_DIR ]
then 
  if ! mkdir -p $RUN_DIR
  then echo "Could not make run dir $RUN_DIR" >&2 ; exit 1
  fi
fi

# Move the run config file in there, if it isn't already
if [ "$RUN_CONFIG_PATH" == "$RUN_DIR/$RUN_CONFIG_NAME" ]
then :
else
  mv $RUN_CONFIG $RUN_DIR
  RUN_CONFIG_PATH=$RUN_DIR/$RUN_CONFIG_NAME
  if [ ! -f "$RUN_CONFIG_PATH" ]
  then echo "Problem placing run config in run dir at $RUN_CONFIG_PATH" >&2; exit 1
  fi
fi
# Reassign RUN_CONFIG to absolute path so we won't lose it when we chdir
RUN_CONFIG=$RUN_CONFIG_PATH
  
# Make and change to the run directory
cd $RUN_DIR
echo "Using run directory $RUN_DIR";

# Get the list of samples on a single line for looping, and save to samples.txt
SAMPLES=`cut -f 1 $RUN_CONFIG | sort | uniq | tee samples.txt`
SAMPLES=`echo $SAMPLES`

# Make the scripts

# Set the run directory
set_step_script_dir $RUN_DIR 

#
# Step 01: BWA
#
init_step bwa
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do make_script_bwa_ubam_rgmerge.sh $GENOME_NAME $DATA_DIR/\$sample/ubam $DATA_DIR/\$sample/align/$GENOME_NAME/WG precopy $ALGO
done > `get_step_cmds_file`
EOF
write_step_submit_script

#
# Step 2: Finish align
#
init_step finishalign
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
for sample in $SAMPLES
do
  dir=$DATA_DIR/\$sample/align/$GENOME_NAME/WG
  bn=\`basename \$dir\`
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$bn-\$sample qc_bam_dir.sh \$dir $DATA_DIR/\$sample/ubam coordinate
  then anyfail=yes
  fi
done
if [ "\$anyfail" == "yes" ]
then 
  echo "There were QA failures"
  echo "Exiting..."
  exit 1
fi
EOF
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
make_script_mmi_from_run_config.sh $RUN_CONFIG $DATA_DIR align/$GENOME_NAME/WG finish > `get_step_cmds_file`
EOF
write_step_submit_script

#
# Step 3: Final QA and stats
# 
init_step finalqa
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# Run final QA step of the Standard Mapping Pipeline

# Run from run dir for simplicity
cd \`dirname \$0\`
# QA:
anyfail=no
while read sample ignored ignored build ignored
do
  dir=$DATA_DIR/\$sample/finish
  bn=\`basename \$dir\`
  if [ "\$build" == "" -o "\$build" == "-" -o "\$build" == "." -o "\$build" == "LATEST" ]
  then buildsuff=
  else buildsuff="-\$build"
  fi
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$bn-\$sample\$buildsuff qc_bam.sh \$dir/\$sample\$buildsuff.bam indexed \$dir/\$sample\$buildsuff.flagstat.txt
  then anyfail=yes
  fi
done < $RUN_CONFIG
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  exit 0
fi
EOF
cat > `get_step_local_work_script` <<EOF
#!/bin/bash
# Compile mapping stats and print them out on stdout

# Stats
echo Mapping stats are written to mapping_stats.txt and are presented here:
( cd $DATA_DIR; compile_mapping_stats.sh -h $SAMPLES ) | tee mapping_stats.txt
awk '{ ndm=\$4; gsub(/,/, "", ndm); print \$1, ndm }' mapping_stats.txt | sample_qc.sh "Mapping Stats" NonDupMapped $CUTOFF
EOF

# Print instructions
echo "Run dir is $RUN_DIR"
