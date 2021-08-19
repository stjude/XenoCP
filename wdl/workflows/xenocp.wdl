## # XenoCP
##
## This WDL workflow runs the XenoCP xenograft cleansing pipeline. 
## The workflow takes an input BAM file and splits it into fastq files for each chromosome. 
## The read pairs are then passed through alignment to generate a BAM file aligned to the host genome.
##
## ## LICENSING
##
## #### MIT License
##
## Copyright 2019 St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, merge,
## publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
## to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

version 1.0

import "https://raw.githubusercontent.com/stjude/xenocp/master/wdl/tools/xenocp.wdl" as xenocp_tools
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/bwa.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/star.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/picard.wdl"

workflow xenocp {
    input {
        File input_bam
        File input_bai
        File reference_tar_gz
        String aligner = "bwa aln"
        String validation_stringency = "SILENT"
        Int n_threads = 1
        Boolean skip_duplicate_marking = false
    }

    parameter_meta {
        input_bam: "BAM file from which to clean contaminate reads"
        reference_tar_gz: "Reference gzipped tar file containing either STAR or BWA indexes depending on which aligner is selected. For BWA, files should be at the root level. For STAR, files should be in a directory that matches the root of the gzipped archive."
        aligner: "Which aligner to use to map reads to the host genome to detect contamination: [bwa aln, bwa mem, star]"
        skip_duplicate_marking: "If true, duplicate marking will be skipped when the cleaned BAMs are merged"
    }
    
    String name = basename(input_bam, ".bam")
    String output_name = name + ".xenocp.bam"

    call xenocp_tools.get_chroms { input: input_bam=input_bam, ncpu=n_threads }
    call xenocp_tools.extract_mismatch as mismatch { input: input_bam=input_bam, input_bai=input_bai }
    scatter (chromosome in get_chroms.chromosomes) {
        call xenocp_tools.extract_by_chrom { input: input_bam=input_bam, input_bai=input_bai, chromosome=chromosome }
    }
    call xenocp_tools.extract_unmapped as unmapped { input: input_bam=input_bam, input_bai=input_bai, ncpu=n_threads }

    Array[File] split_bams = flatten([extract_by_chrom.out_bam, [mismatch.mismatch_bam]])

    scatter (bam in split_bams){
        call xenocp_tools.mapped_fastq { input: input_bam=bam }
    }
    if (aligner == "bwa aln") {
        scatter (fastq in mapped_fastq.fastq){
            call bwa.bwa_aln as bwa_aln_align { input: fastq=fastq, bwadb_tar_gz=reference_tar_gz, ncpu=n_threads }
        }
    }
    if (aligner == "bwa mem") {
        scatter (fastq in mapped_fastq.fastq){
            call bwa.bwa_mem as bwa_mem_align { input: fastq=fastq, bwadb_tar_gz=reference_tar_gz, ncpu=n_threads }
        }
    }
    if (aligner == "star") {
        scatter (fastq in mapped_fastq.fastq){
            call star.alignment as star_align { input: read_one_fastqs=[fastq], stardb_tar_gz=reference_tar_gz, output_prefix=basename(fastq, ".fq.gz"), ncpu=n_threads }
        }
    }
    
    scatter(bam in select_first([bwa_aln_align.bam, bwa_mem_align.bam, star_align.star_bam])){
        call picard.sort as sort { input: bam=bam, sort_order="queryname"}
    }

    scatter (pair in zip(split_bams, sort.sorted_bam)){
        call xenocp_tools.create_contam_list { input: input_bam=pair.left, contam_bam=pair.right, stringency=validation_stringency }
    }

    scatter (pair in zip(split_bams, create_contam_list.contam_list)){
        call xenocp_tools.cleanse { input: input_bam=pair.left, unmap_reads=pair.right, stringency=validation_stringency }
    }

    call xenocp_tools.merge_markdup_index as final_bam { input: input_bams=flatten([cleanse.cleaned_bam, [unmapped.unmapped_bam]]), output_bam=output_name, skip_dup=skip_duplicate_marking }

    call xenocp_tools.qc { input: input_bam=final_bam.final_bam, input_bai=final_bam.final_bai, flagstat=final_bam.flagstat }

    scatter(bam in create_contam_list.output_tie_bam){
        call picard.sort as tie_sort { input: bam=bam, sort_order="coordinate" }
    }

    call xenocp_tools.merge_markdup_index as combine_tie_bam { input: input_bams=tie_sort.sorted_bam, output_bam=name+".tie.bam", skip_dup=true }

    call xenocp_tools.combine_files { input: input_files=create_contam_list.contam_list, output_file=name+".contam.txt" }

    output {
        File bam = final_bam.final_bam
        File bam_index = final_bam.final_bai
        File bam_md5 = final_bam.final_md5
        File flagstat = final_bam.flagstat
        File contam_list = combine_files.final_file
        File tie_bam = combine_tie_bam.final_bam
    }
}
