#!/usr/bin/env bash

### This is a log file to obtain chromosome sizes, FANTOM5 promoters and enhancers, and mask file for enhancer identification. 
### Nov 15 2022

### Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
### Written by:		Shruti Bhagat; bhagat.shruti.6j@kyoto-u.ac.jp; Raku Son; raku.son@riken.jp 
### Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
###					Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 


### Human hg38 chromosome sizes, please visit https://hgdownload.cse.ucsc.edu/goldenpath/ for other species
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.chrom.sizes
sort -k1,1 -k2,2n hg38.chrom.sizes > hg38.chrom.sorted.sizes

### If you use references other than ucsc reference or from other species, chromosome size file can also be produced from the reference fasta by samtools (http://www.htslib.org/)
### An example of creating chromosome size file from GENCODE reference file is shown below
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai > GRCh38.primary_assembly.chrom.sizes

### Human hg38 FANTOM5 promoters, please visit https://fantom.gsc.riken.jp/5/datafiles/ for other CAGE peak files and species
wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
### Human hg38 FANTOM5 enhancers, please visit https://fantom.gsc.riken.jp/5/datafiles/ for other enhancer expression files and species
wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
### Human hg37 FANTOM-NET enhancers, https://fantom.gsc.riken.jp/5/suppl/Hirabayashi_et_al_2019/ for other species
wget https://fantom.gsc.riken.jp/5/suppl/Hirabayashi_et_al_2019/data/Supplementary_Data_1_Human_FANTOM-NET_enhancers.bed.gz

## Unzip
gunzip *.gz

### Get number of promoters and enhancers
for file in *.bed; do
echo $file
wc -l "$file" 
done

### hg38_fair+new_CAGE_peaks_phase1and2.bed
### 209911 
### F5.hg38.enhancers.bed
### 63285
### Supplementary_Data_1_Human_FANTOM-NET_enhancers.bed
### 85786

### Convert bed12 to bed6 by selecting first 6 columns 
for infile in *.bed
 do
 echo ${infile}
   base=$(basename ${infile} .bed)
       awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4, "0", $6 }' ${infile} | sort -k1,1 -k2,2n \
       > ${base}.6.bed
 done

### Edit promoter name to include hg38 coordinates
awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $1":"$2".."$3","$6"|"$4, $5, $6 }' hg38_fair+new_CAGE_peaks_phase1and2.6.bed > hg38_fair+new_CAGE_peaks_phase1and2.edit.bed

### Get number of promoters and enhancers
for file in *.6.bed; do
echo $file
wc -l "$file" 
done

### hg38_fair+new_CAGE_peaks_phase1and2.edit.bed
### 209911 
### F5.hg38.enhancers.6.bed
### 63285 
### Supplementary_Data_1_Human_FANTOM-NET_enhancers.6.bed
### 85786 

### Convert hg37 to hg38 coordinates with default parameters for FANTOM-NET enhancers
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip *.gz
liftOver Supplementary_Data_1_Human_FANTOM-NET_enhancers.6.bed hg19ToHg38.over.chain Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.bed Supplementary_Data_1_Human_FANTOM-NET_enhancers.6.unmapped.bed
wc -l Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.bed
### 85760
wc -l Supplementary_Data_1_Human_FANTOM-NET_enhancers.6.unmapped.bed
### 52
awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $1":"$2".."$3"|hg19::"$4, $5, $6 }' Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.bed > Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.6.bed

### Creating a mask file used for enhancer analysis. Mask transcripts which are not lncRNAs.
### This script uses human GENCODE v39 GTF file. The same can be applied to any other GTF file from the GENCODE database. 

### Remove header from GTF file
awk 'NR > 5 { print }' gencode.v39.primary_assembly.annotation.gtf > Noheader_gencode.v39.primary_assembly.annotation.gtf

### Select all transcripts from gtf file 
awk '$1 ~ /chr/ { print $0 }' Noheader_gencode.v39.primary_assembly.annotation.gtf \
| awk 'BEGIN{OFS="\t"} {if($3 == "transcript") {print $0}}' \
| sort -k1,1 -k4,4n  > Noheader_gencode.v39.chr.transcript.gtf

wc Noheader_gencode.v39.chr.transcript.gtf
### n = 244,939

### Select all non-lncRNA gene-type and transcript-type from gtf file 
awk 'BEGIN{OFS="\t"} {if($14 != "\"lncRNA\";" && $18 != "\"lncRNA\";") {print $0}}' Noheader_gencode.v39.chr.transcript.gtf > Noheader_gencode.v39.chr.gene.transcript.nonlncRNA.gtf

wc Noheader_gencode.v39.chr.gene.transcript.nonlncRNA.gtf
### n = 193,007

### GTF to bed file
awk 'BEGIN{OFS="\t"} {{print $1, $4-1, $5, $1":"$4-1".."$5","$7, "0", $7, $14"-"$18"_"$16} }' Noheader_gencode.v39.chr.gene.transcript.nonlncRNA.gtf > Noheader_gencode.v39.chr.gene.transcript.nonlncRNA.bed

### Extend 5'-end by +/-300bp
awk 'BEGIN{OFS="\t"} {if($6 == "+") {print $1, $2-300, $2+300, $1":"$2-300".."$2+300","$6, $5, $6, $7} else { print $1, $3-300, $3+300, $1":"$3-300".."$3+300","$6, $5, $6, $7} }' Noheader_gencode.v39.chr.gene.transcript.nonlncRNA.bed > Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.bed

### Sort and merge file
sort -k1,1 -k2,2n Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.bed | bedtools merge -s -c 6,7 -o distinct,collapse -delim "|" -i stdin > Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.merged.bed
wc Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.merged.bed
### n = 90,565

awk 'BEGIN{OFS="\t"} {{print $1, $2, $3, $1":"$2".."$3","$4"|"$5, "0", $4} }' Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.merged.bed > Noheader_gencode.v39.chr.gene.transcript.nonlncRNA600.mask.6.bed


