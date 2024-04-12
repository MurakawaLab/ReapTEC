## ReapTEC

This repository provides scripts to run the ReapTEC (read-level pre-filtering and transcribed enhancer call) pipeline. 
By leveraging a unique “cap signature” derived from the 5′-end of a transcript, ReapTEC simultaneously profiles 
gene expression and enhancer activity at nucleotide resolution using 5′-end single-cell RNA-sequencing (5′ scRNA-seq).
The flow of analysis and a brief description of the scripts is provided below. The details are described in the STARsolo_ReapTEC_logfile.

![ReapTEC_pipeline](https://github.com/MurakawaLab/ReapTEC/assets/23185260/126be79d-3788-451a-b9f1-c57fb6d3e9e3)

•	STARsolo_ReapTEC_logfile: This log file provides a step-by-step protocol to run the ReapTEC pipeline. 
All software requirements and other details are provided in this log file. 
The input for ReapTEC is 10x Genomics Chromium paired-end 5′ scRNA-seq fastq files. 
Read 1 needs to be sequenced longer than conventional 26 bp to obtain transcription start sites (TSSs) 
at the 5′-end of transcripts (e.g. we performed 2 × 150 bp paired-end sequencing).

•	Promoter_Enhancer_Download_logFile: This log file provides instructions to obtains reference files, 
FANTOM5 promoters and enhancers, and to create mask file used during enhancer identification.

•	STARsolo_STARindex_GENCODE41_PRI_221224_logFile: This log file provides instructions to create STAR index using GENCODE gene model.

•	STARsolo_SoftclipG_221224.sh: This script extracts reads starting with an unencoded G (at the 5′-end).

•	STARsolo_Cell_barcode_CTSS_bed_20221123.sh: This script identifies 5′-end TSS from STARsolo-aligned umi-tools-deduplicated 5′ scRNA-seq data.

•	STARsolo_Counts_CTSS_bed_CPM_bigwig_240119.sh: This script counts 5′-end of reads that map to known promoters and enhancers. 
This script also provides normalized BigWig files for visualization on IGV, UCSC.

•	fixed_bidir_enhancers_10bp.sh: This script identifies bidirectionally transcribed enhancers. 
This script is adapted from the FANTOM5 consortium (https://github.com/anderssonrobin/enhancers) with minor modifications 
(enhancers/scripts/fixed_bidir_enhancers_10bp.sh).

•	(Optional) TSS_peak_based_analysis_log: This log file provides instructions to perform TSS–based single-cell analysis. 
The 5′ ends of transcripts are counted at the single-cell level for each TSS peak.

