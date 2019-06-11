## Getting started

	git clone https://github.com/adamdingliang/XenoCP.git
	cd XenoCP
	gradle :xenocp:installDist

Download contaminant genomic reference and update `sample_data/input_data/inputs-local.yml` with the path to the reference data.

   
	mkdir results
	cwltool --outdir results cwl/xenocp.cwl sample_data/input_data/inputs_local.yml
	
## Introduction to XenoCP

XenoCP is a cloud-based tool for cleansing mouse reads in xenograft BAMs. XenoCP can be easily incorporated into any workflow as it takes a BAM file
as input and efficiently cleans up the mouse contamination and gives a clean human BAM output that could be used for downstream
genomic analysis. 

St. Jude cloud version

XenoCP workflow:
<!--![Alt text](images/xenocp_workflow2.png) -->
<img src="images/xenocp_workflow2.png" width="500">

## Prerequisites

  * [bwa] =0.7.13
  * [gawk]*
  * [Java SE Development Kit] ~1.8
    * [Gradle] ~5.3
  * [Node.js] ~10.15.3
  * [Picard] =2.6.0
  * [Python] ~3.6.1
    * [cwltool] ~1.0
    * [html5lib] ~1.0.1
  * [samtools]† ~1.9
    * [zlib]

\* XenoCP requires the GNU inplementation of awk.

† XenoCP only supports BAM files. When compiling samtools, CRAM block codecs
(bz2 and lmza) can be disabled. ncurses (for `samtools tview`) can also be
disabled.

[bwa]: https://github.com/lh3/bwa
[gawk]: https://www.gnu.org/software/gawk/
[Java SE Development Kit]: https://www.oracle.com/technetwork/java/javase/overview/index.html
[Gradle]: https://gradle.org/
[Node.js]: https://nodejs.org/en/
[Picard]: https://broadinstitute.github.io/picard/
[Python]: https://www.python.org/
[cwltool]: https://github.com/common-workflow-language/cwltool
[html5lib]: https://github.com/html5lib/html5lib-python
[samtools]: http://www.htslib.org/
[zlib]: https://www.zlib.net/

## Building

Once the prerequisites are satisfied, build XenoCP using Gradle. 

```
$ gradle :xenocp:installDist
```

Add the artifacts under `xenocp/build/install` to your Java CLASSPATH.

## Usage

XenoCP uses [CWL] to describe its workflow.

To run an example workflow, update `sample_data/input_data/inputs_local.yml` with the path to a reference genome.
Then run the following.

```
$ mkdir results
$ cwltool --outdir results cwl/xenocp.cwl sample_data/input_data/inputs_local.yml
```

[CWL]: https://www.commonwl.org/

### Inputs

XenoCP requires three inputs, defined in a YAML file as [CWL inputs]. E.g., `inputs.yml`:

```
bam:
  class: File
  path: sample.bam
bai:
  class: File
  path: sample.bai
ref_db_prefix: /references/ref.fa
```

`bam` is the input sample BAM, `bai` is the bam index for the input sample BAM 
 and `ref_db_prefix`, the basename of the reference assembly that should be cleansed. 
For example, a prefix of `MGSCv37.fa` would assume
the following files in the same directory exist: 
`MGSCv37.fa.amb`, `MGSCv37.fa.ann`, `MGSCv37.fa.bwt`, 
`MGSCv37.fa.pac`, and `MGSCv37.fa.sa`.

Several optional input paramters can be changed in the inputs file.

```
suffix_length: 4
keep_mates_together: true
validation_stringency: SILENT
output_prefix: xenocp-
output_extension: bam
```

### Create Reference Files

Download the FASTA file for your genome assembly and run the following commands to create other files:
```
$ bwa index -p $FASTA $FASTA
```

[CWL inputs]: https://www.commonwl.org/user_guide/02-1st-example/index.html

## Docker

XenoCP provides a [Dockerfile] that builds an image with all the included
dependencies. To use this image, install [Docker] for your platform.

[Docker]: https://www.docker.com/

### Build

In the XenoCP project directory, build the Docker image.

```
$ docker build --tag xenocp .
```

### Run

The Docker image uses `cwl-runner cwl/xenocp.cwl` as its entrypoint.

The image assumes three working directories: `/data` for inputs, `/references` for
reference files, and `/results` for outputs. `/data` and `/references` can be
read-only, where as `/results` needs write access.

The paths given in the input parameters file must be from inside the
container, not the host, e.g.,

```
bam:
  class: File
  path: /data/sample.bam
bai: 
  class: File
  path: /data/sample.bai
ref_db_prefix: /references/ref.fa
```

The following is an example `run` command where files are stored in `test/{data,references}`. Outputs are saved in `test/results`.

This example assumes you are running against Mus Musculus (genome build MGSCv37). Set the path to the folder containing your reference data
and run the following command to produce output from the included sample data. Test output for comparison is located at `sample_data/output_data`.

```
$ mkdir `pwd`/results
$ docker run \
  --mount type=bind,source=`pwd`/sample_data/input_data,target=/data,readonly \
  --mount type=bind,source=/path/to/references,target=/references,readonly \
  --mount type=bind,source=`pwd`/results,target=/results \
  xenocp \
  /data/inputs.yml
```

[Dockerfile]: ./Dockerfile

## St. Jude Cloud

To run XenoCP in St. Jude Cloud, please follow the directions at: https://stjude.github.io/sjcloud-docs/guides/tools/xenocp/

## Availability

[TODO] XenoCP is released under ...

## Seeking help

[TODO]

## Citing XenoCP

[TODO] Publication in prep
