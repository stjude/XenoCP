## # XenoCP
##
## This WDL tool wraps [XenoCP](https://github.com/stjude/XenoCP).
## XenoCP is cleans contamination reads from a BAM file.

version 1.0

task get_chroms {
    input {
        Int ncpu = 1
        File input_bam
        Int memory_gb = 1
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size + 2)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        samtools view -H ~{input_bam} | grep "@SQ" | cut -f 2  | cut -f 2 -d':' > chroms.txt
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File chromosomes_files = "chroms.txt"
        Array[String] chromosomes = read_lines("chroms.txt")
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs TweakSam to extract reads with mate mapped to a different chromosome." 
    }

    parameter_meta {
        input_bam: "Input BAM file from which to extract"
    }
}

task extract_mismatch {
    input {
        Int ncpu = 1
        File input_bam
        File input_bai
        Int memory_gb = 1
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        java.sh org.stjude.compbio.sam.TweakSam -x -O "queryname" -o "other.bam" -i ~{input_bam}

    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File mismatch_bam = "other.bam"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs TweakSam to extract reads with mate mapped to a different chromosome." 
    }

    parameter_meta {
        input_bam: "Input BAM file from which to extract"
    }
}

task extract_by_chrom {
    input {
        Int ncpu = 1
        File input_bam
        File input_bai
        String chromosome
        Int memory_gb = 1
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        java.sh org.stjude.compbio.sam.TweakSam -X -O "queryname" -o "~{chromosome}.bam" -c ~{chromosome} -i ~{input_bam}

    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File out_bam = chromosome + ".bam"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs TweakSam to extract reads mapped to a chromosome." 
    }

    parameter_meta {
        input_bam: "Input BAM file from which to extract"
    }
}

task extract_unmapped {
    input {
        Int ncpu = 1
        File input_bam
        File input_bai
        Int memory_gb = 1
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        sambamba view -h -F "ref_name =~ /\\*/" -f "bam" ~{input_bam} > unmapped.bam

    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File unmapped_bam = "unmapped.bam"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs TweakSam to extract reads with mate mapped to a different chromosome." 
    }

    parameter_meta {
        input_bam: "Input BAM file from which to extract"
    }
}

task mapped_fastq {
    input {
        File input_bam
        String output_fastq = basename(input_bam, ".bam") + ".fastq"
        Int ncpu = 1
        Int memory_gb = 1
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }


    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size + 2)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        view_awk_picard.sh ~{input_bam} ~{output_fastq}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File fastq = output_fastq + ".gz"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool converts a SAM file to FastQ format." 
    }

    parameter_meta {
        input_bam: "Input BAM file from which to extract"
    }
}

task create_contam_list {
    input {
        File input_bam
        String output_contam_list = basename(input_bam, ".bam") + "contam.txt"
        String tie_bam = basename(input_bam, ".bam") + ".tie.bam"
        File contam_bam
        String stringency = "SILENT"
        Int? disk_size_gb
        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 1
    }

    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        java.sh org.stjude.compbio.xenocp.CreateContamLists -i ~{input_bam} -o ~{output_contam_list} -t ~{tie_bam} -V ~{stringency} ~{contam_bam}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File contam_list = output_contam_list
        File output_tie_bam = tie_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool creates a list of contamination reads by comparing host alignments to original alignments." 
    }

    parameter_meta {
        input_bam: "Input BAM file"
    }
}

task cleanse {
    input {
        File input_bam
        String output_bam = basename(input_bam, ".bam") + "-xenocp.bam"
        File unmap_reads
        String sort_order = "coordinate"
        String stringency = "SILENT"
        Int? disk_size_gb
        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 1
    }

    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        java.sh org.stjude.compbio.sam.TweakSam -i ~{input_bam} -o ~{output_bam} -u ~{unmap_reads} -V ~{stringency} -O ~{sort_order}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File cleaned_bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs creates a list of contamination reads by comparing host alignments to original alignments." 
    }

    parameter_meta {
        input_bam: "Input BAM file"
    }
}

task merge_markdup_index {
    input {
        Array[File] input_bams
        Boolean skip_dup = false
        String output_bam = "xenocp.bam"
        Int? disk_size_gb
        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 1
    }

    Float input_bam_size = size(input_bams, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        dup_opt=""
        if [ -n ~{skip_dup} ]
        then
            dup_opt="--no-markdup"
        fi

        merge_markdup_index.sh ${dup_opt} -t ~{ncpu} ~{output_bam} ~{sep=" " input_bams}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
        File final_bam = output_bam
        File final_bai = output_bam + ".bai"
        File final_md5 = output_bam + ".md5"
        File flagstat = basename(output_bam, ".bam") + ".flagstat.txt"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs creates a list of contamination reads by comparing host alignments to original alignments." 
    }

    parameter_meta {
        input_bams: "Input BAM files"
    }
}

task qc {
    input {
        File input_bam
        File input_bai
        String status = "indexed"
        File flagstat
        Int? disk_size_gb
        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 1
    }

    Float input_bam_size = size(input_bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(input_bam_size * 2)])

    command <<<
        set -euo pipefail

        qc_bam.sh ~{input_bam} ~{status} ~{flagstat}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjude/xenocp:latest'
        maxRetries: max_retries
    }

    output {
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs creates a list of contamination reads by comparing host alignments to original alignments." 
    }

    parameter_meta {
        input_bam: "Input BAM file"
    }
}
