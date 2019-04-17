# XenoCP

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

## Usage

XenoCP uses [CWL] to describe its workflow.

```
$ cwltool --outdir tmp/results cwl/xenocp.cwl tmp/inputs.yml
```

[CWL]: https://www.commonwl.org/

### Inputs

XenoCP requires two inputs, defined in a YAML file as [CWL inputs]. E.g., `inputs.yml`:

```
input_bam:
  class: File
  path: sample.bam
ref_db_prefix: /references/ref.fa
```

`input_bam` is the input sample BAM; and `ref_db_prefix`, the basename of the
input reference assembly. For example, a prefix of `MGSCv37.fa` would assume
the following files in the same directory exist: `MGSCv37.fa`,
`MGSCv37.fa.amb`, `MGSCv37.fa.ann`, `MGSCv37.fa.bwt`, `MGSCv37.fa.dict`,
`MGSCv37.fa.fai`, `MGSCv37.fa.pac`, and `MGSCv37.fa.sa`.

Several optional input paramters can be changed in the inputs file.

```
suffix_length: 4
keep_mates_together: true
num_backet: 31
validation_stringency: SILENT
output_prefix: xenocp-
output_extension: bam
sort_order: queryname
```

### Create Reference Files

Download the FASTA file for your genome assembly and run the following commands to create other files:
```
$ samtools faidx $FASTA
$ java picard.cmdline.PicardCommandLine CreateSequenceDictionary REFERENCE=$FASTA OUTPUT=$FASTA.dict
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
input_bam:
  class: File
  path: /data/sample.bam
ref_db_prefix: /references/ref.fa
```

The following is an example `run` command where files are stored in `tmp/test/{data,references}`. Outputs are saved in `tmp/test/results`.

```
$ docker run \
  --mount type=bind,source=$(pwd)/tmp/test/data,target=/data,readonly \
  --mount type=bind,source=$(pwd)/tmp/test/references,target=/references,readonly \
  --mount type=bind,source=$(pwd)/tmp/test/results,target=/results \
  xenocp \
  /data/inputs.yml
```

[Dockerfile]: ./Dockerfile
