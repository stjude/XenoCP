<!-- dx-header -->
# XenoCP: Xenograft sample cleaning pipeline

Analysis of next-generation sequencing (NGS) data from patient-derived
xenograft (PDX) samples is complicated by the presence of reads from
contaminating host cells.  XenoCP cleanses an existing BAM file mapped to the
graft genome by aligning mapped reads to the host genome, identifying reads
that appear to be from host contamination, and marking them as unread in the
BAM file, producing a new cleansed BAM file.

## Host genome selection

By default, the pipeline assumes the host is mouse, and the mouse reference may
be changed in the additional input parameters.

To use another reference genome for host:
1. Prepare a BWA index
2. Tar and gzip the index files into a .tar.gz or .tgz file.
3. Provide it as the input file Custom Host Genome
4. In the additional input parameters, change the Host Genome to Custom.
