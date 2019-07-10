# XenoCP

XenoCP is a tool for cleansing mouse reads in xenograft BAMs.
XenoCP can be easily incorporated into any workflow, as it takes a BAM file
as input and efficiently cleans up the mouse contamination. The output is a clean
human BAM file that could be used for downstream genomic analysis. 

## Getting started

XenoCP can be run in the cloud on DNAnexus at
https://platform.dnanexus.com/app/stjude_xenocp

The easiest way to get XenoCP running locally is using Docker, as `docker build`
creates an image with all of the dependencies:

	git clone https://github.com/stjude/XenoCP.git
	cd XenoCP
	docker build -t xenocp .
	
	wget -r -np -R "index.html*" -nH --cut-dirs=3 http://ftp.stjude.org/pub/software/xenocp/reference/MGSCv37
	
	# Test run on small dataset
	mkdir results
	docker run \
	  --mount type=bind,source=$(pwd)/sample_data/input_data,target=/data,readonly \
	  --mount type=bind,source=$(pwd)/reference,target=/reference,readonly \
	  --mount type=bind,source=$(pwd)/results,target=/results \
	  xenocp \
	  /data/inputs.yml
	
## Introduction to XenoCP

XenoCP takes a BAM with xenograft reads mapped to the graft genome (e.g., human).
It extracts aligned reads and remaps to the host genome (e.g., mouse) to
determine whether the reads are from the host or graft.  The output is a copy of the
original BAM with host reads marked as unmapped.

XenoCP workflow:
<!--![Alt text](images/xenocp_workflow2.png) -->
<img src="images/xenocp_workflow2.png" width="500">

## Prerequisites

  * [bwa] =0.7.13
  * [gawk]*
  * [Java SE Development Kit] ~1.8
    * [Gradle] ~5.3
  * [Node.js] ~10.15.3
  * [Python] ~3.6.1
    * [cwltool] ~1.0
    * [html5lib] ~1.0.1
  * [samtools]† ~1.9
    * [zlib]
  * [sambamba] ~0.7.0

\* XenoCP requires the GNU inplementation of awk.

† XenoCP only supports BAM files. When compiling samtools, CRAM block codecs
(bz2 and lmza) can be disabled. ncurses (for `samtools tview`) can also be
disabled.

[bwa]: https://github.com/lh3/bwa
[gawk]: https://www.gnu.org/software/gawk/
[Java SE Development Kit]: https://www.oracle.com/technetwork/java/javase/overview/index.html
[Gradle]: https://gradle.org/
[Node.js]: https://nodejs.org/en/
[Python]: https://www.python.org/
[cwltool]: https://github.com/common-workflow-language/cwltool
[html5lib]: https://github.com/html5lib/html5lib-python
[samtools]: http://www.htslib.org/
[zlib]: https://www.zlib.net/
[sambamba]: http://lomereiter.github.io/sambamba/



## Local usage


### Obtain XenoCP 

Clone XenoCP from GitHub: 
```
git clone https://github.com/stjude/XenoCP.git
```

### Build XenoCP

Once the prerequisites are satisfied, build XenoCP using Gradle. 

```
$ gradle installDist
```

Add the artifacts under `build/install/xenocp/lib` to your Java `CLASSPATH`.
Add the artifacts under `build/install/xenocp/bin` to your `PATH`.

### Inputs

XenoCP requires two inputs, defined in a YAML file as [CWL inputs]. E.g., `inputs.yml`:

```
bam:
  class: File
  path: sample.bam
ref_db_prefix: /references/ref.fa
```

`bam` is the input sample BAM 
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

### Download MGSCv37 reference files

Reference files are provided for version MGSCv37 of mouse and are available from http://ftp.stjude.org/pub/software/xenocp/reference/MGSCv37

### Run

XenoCP uses [CWL] to describe its workflow.

To run an example workflow, update `sample_data/input_data/inputs_local.yml` with the path to a reference genome.
Then run the following.

```
$ mkdir results
$ cwltool --outdir results cwl/xenocp.cwl sample_data/input_data/inputs_local.yml
```

[CWL]: https://www.commonwl.org/

## Docker

XenoCP provides a [Dockerfile] that builds an image with all the included
dependencies. To use this image, install [Docker] for your platform.

[Docker]: https://www.docker.com/

### Build Docker image

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
ref_db_prefix: /reference/ref.fa
```

The following is an example `run` command where files are stored in `test/{data,reference}`. Outputs are saved in `test/results`.

This example assumes you are running against Mus musculus (genome build MGSCv37). Set the path to the folder containing your reference data
and run the following command to produce output from the included sample data. Test output for comparison is located at `sample_data/output_data`.

```
$ mkdir $(pwd)/results
$ docker run \
  --mount type=bind,source=$(pwd)/sample_data/input_data,target=/data,readonly \
  --mount type=bind,source=/path/to/reference,target=/reference,readonly \
  --mount type=bind,source=$(pwd)/results,target=/results \
  xenocp \
  /data/inputs.yml
```

[Dockerfile]: ./Dockerfile

## Evaluate test data results

If you have [bcftools] and a [GRCh37-lite] reference file, the following will show two variants in the input file. 
The variant on chromosome 1 is an artifact of mouse reads. The variant on chromosome 9 is a variant in the graft genome. 

```
$ bcftools mpileup -R sample_data/output_data/regions.bed -f ref/GRCh37-lite/GRCh37-lite.fa sample_data/input_data/SJRB001_X.subset.bam | bcftools call -m - | tail -n 3
```

Output: 
```
[mpileup] 1 samples in 1 input files
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_LC-SJRB001-X-SJ39.3-8L
1	156044156	.	C	T	67	.	DP=49;VDB=0.525878;SGB=-0.670168;RPB=0.999118;MQB=0.00218214;MQSB=0.948436;BQB=0.743365;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=15,19,2,8;MQ=52	GT:PL	0/1:102,0,255
9	19451994	.	G	A	182	.	DP=26;VDB=0.130558;SGB=-0.680642;RPB=0.887078;MQB=0.948139;MQSB=0.955682;BQB=0.053431;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=5,8,6,6;MQ=58	GT:PL	0/1:215,0,255
```

After running XenoCP, the host genome variant is removed, as the supporting reads will be unmapped. The following command demonstrates the removal of the variant on chromosome 1 in the output of the sample data.


```
$ bcftools mpileup -R sample_data/output_data/regions.bed -f ref/GRCh37-lite/GRCh37-lite.fa sample_data/output_data/SJRB001_X.subset.xenocp.bam | bcftools call -m - | tail -n 3
```

Output: 
```
[mpileup] 1 samples in 1 input files
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_LC-SJRB001-X-SJ39.3-8L
1	156044156	.	C	.	285	.	DP=41;MQSB=0.633762;MQ0F=0;AN=2;DP4=16,20,0,0;MQ=57	GT	0/0
9	19451994	.	G	A	191	.	DP=27;VDB=0.198993;SGB=-0.683931;RPB=0.729125;MQB=0.945959;MQSB=0.960078;BQB=0.0425475;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=5,8,6,7;MQ=58	GT:PL	0/1:224,0,255
```


[bcftools]: https://samtools.github.io/bcftools/bcftools.html
[GRCh37-lite]: ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz


## St. Jude Cloud

To run XenoCP in St. Jude Cloud, please follow the directions at https://www.stjude.cloud/docs/guides/tools/xenocp/

## Availability

Copyright 2019 St. Jude Children's Research Hospital

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Seeking help

For questions and bug reports, please open an issue on the [GitHub] project page.

[GitHub]: https://github.com/stjude/XenoCP/issues

## Citing XenoCP

[TODO] Publication in prep

## Common Issues

Sambamba uses a large number of temporary files while merging the final bam file. Depending on your system, the default open file limit may be too low. You can check the limit with `ulimit -n` and set the limit higher with `ulimit -Sn <value>`.
