#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 November 23rd 2022
 This script can be used to gain 5'-end bed files from STARsolo-aligned umi-tools-deduplicated 5'scRNA-seq data.
 All input BAM files (only Read1 if you have used paired-end data) must be in the working directory and have the suffix ".bam".
  
usage: $0  

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Raku Son; raku.son@riken.jp 
                Akiko Oguchi; akiko.oguchi@riken.jp   
Reviewed by:   	Sho Sekito; sekito.shou.25u@st.kyoto-u.ac.jp
                Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
                Yasuhiro Murakawa; yasuhiro.murakawa@riken.jp, murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}

for infile in *bam; 
 do base=$(basename ${infile} .bam)
 bamToBed -i ${infile} \
  | paste - ${base}_cell_barcode_tmp.txt \
  | awk 'BEGIN {OFS="\t"} { 
    if($6=="+"){print $1, $2, $2+1, ".", $7, $6} 
    else {print $1, $3-1, $3, ".", $7, $6} 
    }' \
  | sort -k1,1 -k2,2n -k6,6 \
  | bedtools groupby -g 1,2,3,4,5,6 -c 1 -o count \
  | awk 'BEGIN{OFS="\t"}{if ($1 ~ /chr/) print $1, $2, $3, $5, $7, $6}' \
  > ${base}.CTSS.fwd.rev.cell.barcode.bed
done
