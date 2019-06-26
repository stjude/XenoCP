#!/bin/bash
# stjude_xenocp
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.
set -e -o pipefail
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

main() {
    echo "Value of input_bam: '$input_bam'"
    echo "Value of input_bai: '$input_bai'"
    echo "Value of ref_name: '$ref_name'"
    echo "Value of ref_bwa_index_tgz: '$ref_bwa_index_tgz'"
    echo "Value of output_prefix: '$output_prefix'"
    echo "Value of output_extension: '$output_extension'"
    
    # This tool uses three directories: data, reference, output, and each has
    # to exist locally and also be bound to a path within the container.  All
    # 3 x 2=6 paths are set here
    # Note: local_output_dir has to be $HOME/out for dx-upload-all-outputs
    # Note: container_output_dir has to be /results because the docker
    # entry point specifies that as the output directory for the cwl command
    local_data_dir=$HOME/data
    local_reference_dir=$HOME/reference
    local_output_dir=$HOME/out
    container_data_dir=/data
    container_reference_dir=/reference
    container_output_dir=/results
    
    echo ""
    echo "=== Setup ==="
    echo "  [*] Downloading input files ..."
    mkdir $local_data_dir
    dx download "$input_bam" -o $local_data_dir/input.bam
    dx download "$input_bai" -o $local_data_dir/input.bam.bai
    
    # We allow the user to choose common host genomes from a dropdown
    # (ref_name), with one option being Custom, which needs to be accompanied
    # by a tar-gzipped bwa index being specified as an input
    # (ref_bwa_index_tgz).  Here, we pull the bwa index from the appropriate
    # place and determine the file basename.  If neither or both of the two
    # options are specified, we error out, as the index cannot be
    # unambiguously determined in those cases.
    echo "  [*] Downloading reference files ..."
    mkdir $local_reference_dir
    if [ "$ref_name" != "" -a ${ref_name:0:6} != "Custom" ]
    then
      if [ "$ref_bwa_index_tgz" != "" ]
      then
        echo "Could not determine which host genome to use." >&2
        echo "A custom host genome was provided but Host Genome was not set to custom." >&2
        exit 1
      fi
      dx download -o $local_reference_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/Mus_musculus/$ref_name/BWA
      mv $local_reference_dir/BWA/* $local_reference_dir/
    else
      if [ "$ref_bwa_index_tgz" == "" ]
      then
        echo "Could not determine which host genome to use." >&2
        echo "A custom host genome was not provided but Host Genome was set to custom." >&2
        exit 1
      fi
      dx download -o $local_reference_dir $ref_bwa_index_tgz
      (
        cd $local_reference_dir
        tar xzf *
      )
    fi
    reference_prefix=`cd $local_reference_dir; echo *.bwt | sed 's/.bwt$//'`
    if ! ls $local_reference_dir/$reference_prefix*
    then echo "Problem obtaining bwa index" >&2; exit 1
    fi

    echo "  [*] Printing instance information ..."
    cat /etc/lsb-release

    echo "  [*] Creating directories and inputs.yml ..."
    mkdir -p $local_output_dir/output_bam
    mkdir -p $local_output_dir/flagstat
    mkdir -p $local_output_dir/contam_list
    mkdir -p $local_output_dir/output_tie_bam
    cat > $local_data_dir/inputs.yml <<EOF
bam:
  class: File
  path: $container_data_dir/input.bam
ref_db_prefix: $container_reference_dir/$reference_prefix
output_prefix: $output_prefix
output_extension: $output_extension
EOF

    echo "  [*] Loading container image ..."
    image_tarfile_path=/stjude/xenocp-docker.tar
    if [ -e $image_tarfile_path.gz ]
    then gunzip $image_tarfile_path.gz
    fi
    docker load -i $image_tarfile_path
    
    echo "=== Execution ==="
    
    # Don't make assumptions about the tag that was used when the image was
    # built, other than that it should be "xenocp".  Since this is run in a
    # clean environment, and we only did a single docker load, the method
    # below should be a safe way to determine the image.
    echo "  [*] Determine image ID ..."
    image_id=`docker images -q xenocp | head -n 1`
    
    echo "  [*] Run container ..."
    docker run \
      --mount type=bind,source=$local_data_dir,target=$container_data_dir,readonly \
      --mount type=bind,source=$local_reference_dir,target=$container_reference_dir,readonly \
      --mount type=bind,source=$local_output_dir,target=$container_output_dir \
      $image_id \
      $container_data_dir/inputs.yml

    echo "=== Wrap Up ==="
    echo "  [*] Uploading outputs ..."
    mv -v $local_output_dir/*.xenocp.bam* $local_output_dir/output_bam
    mv -v $local_output_dir/*.xenocp.flagstat.txt $local_output_dir/flagstat
    mv -v $local_output_dir/*.contam.txt $local_output_dir/contam_list
    mv -v $local_output_dir/*.tie.bam $local_output_dir/output_tie_bam
    ls -l $local_output_dir/*
    dx-upload-all-outputs
    echo "All done!"
}
