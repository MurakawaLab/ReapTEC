## ReapTEC

This repository provides scripts to run the ReapTEC (read-level pre-filtering and transcribed enhancer call) pipeline. 
By leveraging a unique “cap signature” derived from the 5′ end of a transcript, ReapTEC simultaneously profiles 
gene expression and enhancer activity at nucleotide resolution using 5′-end single-cell RNA-sequencing (5′ scRNA-seq).
The flow of analysis and a brief description of the scripts are provided below. The details are described in the STARsolo_ReapTEC_logfile.
All scripts are in the bin directory, please check the individual scripts for detailed usage. A brief description is provided below:


<img width="1539" alt="Screenshot 2024-04-12 at 12 08 03" src="https://github.com/MurakawaLab/ReapTEC/assets/23185260/efd00e63-e106-4d2a-b249-832400d132a2">


•	STARsolo_ReapTEC_logfile: This log file provides a step-by-step protocol to run the ReapTEC pipeline. 
All software requirements and other details are provided in this log file. 
The input for ReapTEC is 10x Genomics Chromium paired-end 5′ scRNA-seq fastq files. 
Read 1 needs to be sequenced longer than conventional 26 bp to obtain transcription start sites (TSSs) 
at the 5′ end of transcripts (e.g. we performed 2 × 150 bp paired-end sequencing).

•	Promoter_Enhancer_Download_logFile: This log file provides instructions to obtain reference files, 
FANTOM5 promoters and enhancers, and to create mask file used during enhancer identification.

•	STARsolo_STARindex_GENCODE41_PRI_221224_logFile: This log file provides instructions to create STAR index using GENCODE gene model.

•	STARsolo_SoftclipG_221224.sh: This script extracts reads starting with an unencoded G (at the 5′ end).

•	STARsolo_Cell_barcode_CTSS_bed_20221123.sh: This script identifies 5′-end TSS from STARsolo-aligned umi-tools-deduplicated 5′ scRNA-seq data.

•	STARsolo_Counts_CTSS_bed_CPM_bigwig_240119.sh: This script counts 5′ end of reads that map to known promoters and enhancers. 
This script also provides normalized BigWig files for visualization on IGV, UCSC.

•	fixed_bidir_enhancers_10bp.sh: This script identifies bidirectionally transcribed enhancers. 
This script is adapted from the FANTOM5 consortium (https://github.com/anderssonrobin/enhancers) with minor modifications 
(enhancers/scripts/fixed_bidir_enhancers_10bp.sh).

•	(Optional) TSS_peak_based_analysis_log: This log file provides instructions to perform TSS-based single-cell analysis. 
The 5′ ends of transcripts are counted at the single-cell level for each TSS peak.

Related data is deposited here https://doi.org/10.5061/dryad.gtht76hv9 and described in the README_Dryad_for_github_April15_2024 file. 

#### Note: The ReapTEC scripts are applicable to data from Next GEM versions. Our team is currently testing whether ReapTEC can be applied to GEM-X data as well.

#### Note: 
When using ReapTEC for 5′ GEM-X, please change the following parts as shown below.

**For STAR mapping:**
--soloUMIstart 17 --soloUMIlen 12   #← change "12"
--clip5pNbases 41 0                #← change "41"

**For softclipG:**
BASE = substr($10, 42, 1);        # ← change "42"
if ($6 ~ /^42S[0-9]/ && BASE == "G") {print $0} \　#← change "42"
ALT = substr($10, length($10)-41, 1); #← change "41"
if ($6 ~ /[0-9]M42S$/ && ALT == "C") {print $0} \ #← change "42"
