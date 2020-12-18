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

import "../../tools/xenocp.wdl" as xenocp_tools
import "../../tools/bwa.wdl"

workflow xenocp {
    input {
        File input_bam
        File input_bai
        File bwadb_tar_gz
        String aligner = "bwa aln"
        Int suffix_length = 4
        Boolean keep_mates_together = true
        String validation_stringency = "SILENT"
        String output_prefix = "xenocp-"
        String output_extension = "bam"
        Int n_threads = 1
    }

    parameter_meta {

    }

    String name = basename(input_bam, ".bam") + ".xenocp.bam"

    call xenocp_tools.get_chroms { input: input_bam=input_bam }
    call xenocp_tools.extract_mismatch as mismatch { input: input_bam=input_bam, input_bai=input_bai }
    scatter (chromosome in get_chroms.chromosomes) {
        call xenocp_tools.extract_by_chrom { input: input_bam=input_bam, input_bai=input_bai, chromosome=chromosome }
    }
    call xenocp_tools.extract_unmapped as unmapped { input: input_bam=input_bam, input_bai=input_bai }

    Array[File] split_bams = flatten([extract_by_chrom.out_bam, [mismatch.mismatch_bam]])

    scatter (bam in split_bams){
        call xenocp_tools.mapped_fastq { input: input_bam=bam }
    }

    scatter (fastq in mapped_fastq.fastq){
        call bwa.bwa_aln as align { input: fastq=fastq, bwadb_tar_gz=bwadb_tar_gz }
    }

    scatter (pair in zip(split_bams, align.bam)){
        call xenocp_tools.create_contam_list { input: input_bam=pair.left, contam_bam=pair.right }
    }

    scatter (pair in zip(split_bams, create_contam_list.contam_list)){
        call xenocp_tools.cleanse { input: input_bam=pair.left, unmap_reads=pair.right }
    }

    call xenocp_tools.merge_markdup_index { input: input_bams=flatten([cleanse.cleaned_bam, [unmapped.unmapped_bam]]), output_bam=name }

    call xenocp_tools.qc { input: input_bam=merge_markdup_index.final_bam, input_bai=merge_markdup_index.final_bai, flagstat=merge_markdup_index.flagstat }

    output {
        File bam = merge_markdup_index.final_bam
        File bam_index = merge_markdup_index.final_bai
        File bam_md5 = merge_markdup_index.final_md5
        File flagstat = merge_markdup_index.flagstat
        Array[File] contam_list = create_contam_list.contam_list
        Array[File] tie_bam = create_contam_list.output_tie_bam
    }
}