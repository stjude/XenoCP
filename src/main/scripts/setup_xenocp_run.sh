#!/bin/bash
# Sets up a run of XenoCP
#
# Directory structure:
# <data directory>/<sample>/xenocp
# +-- orig
#     +-- <sample>.bam (input)
#     +-- <sample>.flagstat.txt (input)
# +-- orig-split-unsorted
#     +-- xenocp-<nnnn>.bam
# +-- orig-split
#     +-- xenocp-<nnnn>.bam
#     +-- xenocp-<nnnn>.flagstat.txt
# +-- fastq
#     +-- xenocp-<nnnn>.bam
# +-- align
#     +-- <genome>
#         +-- xenocp-<nnnn>.bam
# +-- contam
#     +-- xenocp-<nnnn>.contam.txt
# +-- tie
#     +-- xenocp-<nnnn>.tie.bam
# +-- cleansed-split
#     +-- xenocp-<nnnn>.bam
# +-- cleansed (default output)
#     +-- <sample>[-<qualifier>].bam
# <output directory> (output, if specified)
# +-- <sample>
#     +-- <target>
#         +-- <sample>[-<qualifier>].bam
#
# The run config file is tab-delimited text with the following columns:
# 1. Full path to the BAM file
# 2. Sample name
# 3. Target name
# 4. Qualifier (or blank/omitted for unqualified)
#
# Parameters:
# $1 = contaminant genome(s); if multiple, use a comma-delimited list, e.g. X,Y
# $2 = run configuration file
# $3 = run directory
# $4 = data directory
# $5 = (optional) output data directory, if different than main data dir

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
GENOMES="$1"
RUN_CONFIG=$2
RUN_DIR=$3
DATA_DIR=$4
OUT_DIR=$5

# Get pipeline config
. import_config.sh app xenocp

# Validate run config
if [ ! -f $RUN_CONFIG ]
then echo "Run config not found: $RUN_CONFIG" >&2; exit 1
fi
if [ `cut -f 2 $RUN_CONFIG | sort | uniq | wc -l` != `cut -f 2-4 $RUN_CONFIG | sort | uniq | wc -l` ]
then echo "Only one target/qualifier per sample is supported in each run" >&2; exit 1
fi

# Convert genomes to space-delimited
GENOMES=`echo $GENOMES | tr , ' '`

# Create directories if missing
mkdir -p $DATA_DIR
if [ "$OUT_DIR" != "" ]; then mkdir -p $OUT_DIR; fi

# Import step library
. steplib.sh

# Set the run directory
set_step_script_dir $RUN_DIR

# Get the list of samples on a single line for looping, and save to samples.txt
mkdir -p $RUN_DIR
SAMPLES=`cut -f 2 $RUN_CONFIG | sort | uniq | tee $RUN_DIR/samples.txt`
SAMPLES=`echo $SAMPLES`

#
# Step 1: Organize and get flagstats
#
init_step init
# First, make symlinks in the canonical location, if necessary, as a local step
cat > `get_step_local_work_script` <<EOF
#!/bin/bash
exitcode=0
while read bam sample target qualifier
do
	dir=$DATA_DIR/\$sample
	mkdir -p \$dir
	if [ "\$bam" != "" -a ! -e \$dir/\$sample.bam ]
	then ln -s \$bam \$dir/\$sample.bam
	fi
	if [ ! -s \$dir/\$sample.bam ]
	then echo "Input not found at \$dir/\$sample.bam (bam in config was \$bam)" >&2; exitcode=1
	fi
	fs=\`dirname \$bam\`/\`basename \$bam .bam\`.flagstat.txt
	if [ -s "\$fs" -a ! -s "\$dir/\$sample.flagstat.txt" ]
	then ln -s \$fs \$dir/\$sample.flagstat.txt
	fi
done < $RUN_CONFIG
exit \$exitcode
EOF
# Now, create missing flagstats
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  dir=$DATA_DIR/\$sample
  if [ ! -e \$dir/\$sample.flagstat.txt ]
  then echo samtools flagstat \$dir/\$sample.bam \> \$dir/\$sample.flagstat.txt
  fi
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 2: Split input file
#
init_step split
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  input=\$rootdir/\$sample.bam
  size=\`stat -L -c '%s' \$input\`
  rgs=\`samtools view -H \$input | grep -c @RG\`
  b=\`expr \$size / \$rgs / 300000000 + 1\`
  mkdir -p \$rootdir/orig-split-unsorted
  echo java.sh org.stjude.compbio.sam.SplitSam -i \$input -a 4 -m -b \$b -V SILENT \$rootdir/orig-split-unsorted/xenocp- bam
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 3: Sorting split files
#
init_step orig-sort
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  for input in \$rootdir/orig-split-unsorted/*.bam
  do echo sort_flagstat.sh \$input queryname \$rootdir/orig-split/\`basename \$input\`
  done
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 4: Extract FASTQ for mapped reads
#
init_step mapped-fastq
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$sample qc_bam_split.sh \$rootdir/orig-split \$rootdir/\$sample.flagstat.txt
  then anyfail=yes
  fi
done
if [ "\$anyfail" == "yes" ]
then 
  echo "There were QC failures"
  echo "Exiting..."
  exit 1
fi
EOF
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  make_script_xeno_fastq.sh \$rootdir/orig-split \$rootdir/fastq
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 5: Run BWA against contaminant genomes
#
init_step mapping
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for db in $GENOMES
do
  for sample in $SAMPLES
  do
    rootdir=$DATA_DIR/\$sample
    make_script_bwa_alignse_onlymapped.sh \$db \$rootdir/fastq \$rootdir/align/\$db
  done
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 6: Find contamination
#
init_step contamination
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  make_script_find_contam.sh \$rootdir/orig-split \$rootdir/contam \$rootdir/tie \$rootdir/align/*
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 7: Cleanse BAMs
#
init_step cleanse
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
for sample in $SAMPLES
do
  rootdir=$DATA_DIR/\$sample
  make_script_set_reads_unmapped.sh \$rootdir/orig-split \$rootdir/contam \$rootdir/cleansed-split
done > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 8: finish (merge, markdup, index)
#
init_step finish
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
while read bam sample target qualifier
do
	if [ "$OUT_DIR" == "" ]
	then dir=$DATA_DIR/\$sample/xenocp/cleansed
	else dir=$OUT_DIR/\$sample.\$qualifier/\$target
	fi
	mkdir -p \$dir
	if [ "\$qualifier" == "" ]; then dotq= ; else dotq=".\$qualifier"; fi
	echo merge_markdup_index.sh \$dir/\$sample\$dotq.bam $DATA_DIR/\$sample/cleansed-split/*.bam
done < $RUN_CONFIG > `get_step_cmds_file`
EOF
write_step_submit_script


#
# Step 9: Final QC
# 
init_step finalqc
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
anyfail=no
while read bam sample target qualifier
do
	if [ "$OUT_DIR" == "" ]
	then dir=$DATA_DIR/\$sample/xenocp/cleansed
	else dir=$OUT_DIR/\$sample.\$qualifier/\$target
	fi
	if [ "\$qualifier" == "" ]; then dotq= ; else dotq=".\$qualifier"; fi
	if ! qcquiet.sh --exit-on-warn 0 `get_step_failed_qc_dir`/\$sample qc_bam.sh \$dir/\$sample\$dotq.bam indexed \$dir/\$sample\$dotq.flagstat.txt
	then anyfail=yes
	fi
done < $RUN_CONFIG
if [ "\$anyfail" == "yes" ]
then 
  echo There were QC failures
  echo "Exiting..."
  exit 1
fi
EOF


# Print run dir
echo "Run dir is $RUN_DIR"
