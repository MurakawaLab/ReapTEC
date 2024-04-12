#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 December 24th 2022
 This script filters 5' scRNA-seq bam files (mapped by STARsolo) by only retaining reads starting with a soft-clipped G.
 
 This script uses Samtools-1.14 or higher (http://www.htslib.org/).
 
 The requirement is as follows:
 -t Thread count.
 
usage: $0 -t integer 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Raku Son; raku.son@riken.jp 
                Akiko Oguchi; akiko.oguchi@riken.jp   
Reviewed by:   	Sho Sekito; sekito.shou.25u@st.kyoto-u.ac.jp
                Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
                Yasuhiro Murakawa; yasuhiro.murakawa@riken.jp, murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}

while getopts t: opt
do
  case ${opt} in
  t) IntegerT=${OPTARG};;
  *) usage;;
  esac
done

if [ "${IntegerT}" = "" ]; then usage; fi

BASE=G
ALT=C

for infile in *bam
 do
  base=$(basename ${infile} .bam)
  samtools view -@ ${IntegerT} -H ${infile} > ${base}_header.sam

### Gain reads with one-nucleotide soft-clipped on the forward strand
  samtools view -@ ${IntegerT} -F 16 ${infile} | awk -F '\t' '
  BEGIN {OFS="\t"} {
  BASE = substr($10, 40, 1);
  if ($6 ~ /^40S[0-9]/ && BASE == "G") {print $0} \
  }
 ' \
  > ${base}_SoftclipG_F.sam

 ### Gain reads with one-nucleotide soft-clipped on the reverse strand
  samtools view -@ ${IntegerT} -f 16 ${infile} | awk -F '\t' '
  BEGIN {OFS="\t"} {
  ALT = substr($10, length($10)-39, 1);
    if ($6 ~ /[0-9]M40S$/ && ALT == "C") {print $0} \
    }
   ' \
  >  ${base}_SoftclipG_R.sam

### Combine the header, F.sam and R.sam
cat \
  ${base}_header.sam \
  ${base}_SoftclipG_F.sam \
  ${base}_SoftclipG_R.sam \
 | samtools sort -@ ${IntegerT} -O bam -o SoftclipG_${base}.bam

### Delete unnecessary files
rm ${base}*.sam
done