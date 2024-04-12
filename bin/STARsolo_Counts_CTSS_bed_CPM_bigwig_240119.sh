#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 Jan 19th 2024
 This script can be used to gain promoter and enhancer counts and depth normalized BigWig files for visualization on IGV, UCSC. 
 
 This script uses BEDTools and can be installed from here: https://github.com/arq5x/bedtools2.
 bedGraphToBigWig can be downloaded from here: http://hgdownload.soe.ucsc.edu/admin/exe/.
 Please edit your PATH environmental variable in bash_profile so that the above tools can be executable from any directory.
 
 And /path/to/files/ need be provided for the requirements given below:
 -G Sorted chrom sizes file (Human hg38 can be found here: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/ The UCSC goldenpath also contains referance files for all other species).
 -p A bed 6 file containing promoters (For example, FANTOM5 hg38 data can be found here: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz please select the fisrt 6 columns using awk).
 -e A bed 6 file containing enhancers (For example, FANTOM-NET hg37 data can be found here: please use liftOver to convert it to hg38: https://fantom.gsc.riken.jp/5/suppl/Hirabayashi_et_al_2019/data/Supplementary_Data_1_Human_FANTOM-NET_enhancers.bed.gz or FANTOM5 hg38 can be found here: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz please select the fisrt 6 columns using awk).
 All input BED files (only containing 1-base long regions) must be in the working directory and have the suffix ".CTSS.fwd.rev.cell.barcode.bed".
 
 The output files are as follows:
 PREFIX.promoter.fwd.rev.txt
 PREFIX.enhancer.fwd.txt
 PREFIX.enhancer.rev.txt
 PREFIX.CTSS.CPM.fwd.bw
 PREFIX.CTSS.CPM.rev.bw
 
usage: $0 -G genome -p promoter.bed -e enhancer.bed 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp; Raku Son; raku.son@riken.jp 
Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
		Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
		Yasuhiro Murakawa; yasuhiro.murakawa@riken.jp, murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}


while getopts G:p:e: opt
do
  case ${opt} in
  G) genome=${OPTARG};;
  p) promoter=${OPTARG};;
  e) enhancer=${OPTARG};;
  *) usage;;
  esac
done


if [ "${genome}" = "" ]; then usage; fi
if [ "${promoter}" = "" ]; then usage; fi
if [ "${enhancer}" = "" ]; then usage; fi


### Ensure that column5 of ${promoter} and ${enhancer} to be 0
awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, "0", $6}' ${promoter} > Promoter.tmp.bed
awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, "0", $6}' ${enhancer} > Enhancer.tmp.bed

### Counting TSS for promoter
echo "Counting TSS for promoter and enhancer"
for infile in *.CTSS.fwd.rev.cell.barcode.bed
do
   base=$(basename ${infile} .CTSS.fwd.rev.cell.barcode.bed)
   echo ${infile}
   cat ${infile} Promoter.tmp.bed \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge -d -1 -s -c 5,6 -o sum,distinct -i stdin \
  | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, ".", $4, $5}' - \
  | bedtools intersect -wa -wb -s -a Promoter.tmp.bed -b stdin \
  | cut -f 4,11 | sort \
  > ${base}.promoter.fwd.rev.txt

### Counting TSS for enhancer on the forward strand
    awk '{if($6 == "+"){print}}' ${infile} \
  | cat - Enhancer.tmp.bed \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge -d -1 -c 5 -o sum -i stdin \
  | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, ".", $4, "+"}' - \
  | bedtools intersect -wa -wb -a Enhancer.tmp.bed -b stdin \
  | cut -f 4,11 | sort \
  > ${base}.enhancer.fwd.txt

### Counting TSS for enhancer on the reverse strand
    awk '{if($6 == "-"){print}}' ${infile} \
  | cat - Enhancer.tmp.bed \
  | sort -k 1,1 -k 2,2n \
  | bedtools merge -d -1 -c 5 -o sum -i stdin \
  | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, ".", $4, "-"}' - \
  | bedtools intersect -wa -wb -a Enhancer.tmp.bed -b stdin \
  | cut -f 4,11 | sort \
  > ${base}.enhancer.rev.txt
 done


### CPM transform CTSS bed files and convert to bigWig for visualization
echo "Producing CPM normalized bigwig for visualization"
for infile in *.CTSS.fwd.rev.cell.barcode.bed
 do
   echo ${infile}
   base=$(basename ${infile} .CTSS.fwd.rev.cell.barcode.bed)   
   ### Sum CTSS 
   sum=$(awk 'BEGIN{sum=0}{sum=sum+$5}END{print sum}' ${infile} )
   ### Collapse CTSS bed files to avoid overlapping regions in bedGraph
   awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $5, $6}' ${infile} \
   | bedtools groupby -g 1,2,3,4,6 -c 5 -o sum \
   | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $6, $5}' > ${base}.collapsed.bed

   ### CPM transform forward
   awk --assign sum=${sum} 'BEGIN{OFS="\t"} { if($6=="+") { printf("%s\t%i\t%i\t%1.2f\n", $1,$2,$3, 1e6 * $5 / sum) } }' ${base}.collapsed.bed \
   | sort -k1,1 -k2,2n > ${base}.CTSS.CPM.fwd.bg
   
   ### CPM transform reverse
   awk --assign sum=${sum} 'BEGIN{OFS="\t"} { if($6=="-") { printf("%s\t%i\t%i\t%1.2f\n", $1,$2,$3, 1e6 * $5 / sum) } }'  ${base}.collapsed.bed \
   | sort -k1,1 -k2,2n > ${base}.CTSS.CPM.rev.bg
   
   ### Convert to bigWig forward
   bedGraphToBigWig ${base}.CTSS.CPM.fwd.bg ${genome} ${base}.CTSS.CPM.fwd.bw
   
   ### Convert to bigWig reverse
   bedGraphToBigWig ${base}.CTSS.CPM.rev.bg ${genome} ${base}.CTSS.CPM.rev.bw
 done

### Delete unnecessary files
rm Promoter.tmp.bed
rm Enhancer.tmp.bed
rm *.collapsed.bed
rm *.CTSS.CPM.fwd.bg
rm *.CTSS.CPM.rev.bg

### Arrange directories
mkdir Counts_all_cells
mkdir CPMbw_all_cells
mv *.promoter.fwd.rev.txt Counts_all_cells
mv *.enhancer.fwd.txt Counts_all_cells
mv *.enhancer.rev.txt Counts_all_cells
mv *.CPM.* CPMbw_all_cells
