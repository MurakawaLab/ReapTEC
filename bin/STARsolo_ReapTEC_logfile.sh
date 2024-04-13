### This is a log file to obtain bidirectionally transcribed enhancers from 10x Genomics Chromium paired-end 5' scRNA-seq data (Read 1 needs to be sequenced longer than 26bp to obtain 5'-end of transcripts. e.g. we performed 2 × 150 bp paired-end sequencing.)
### April 29th 2023
### Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
### Written by: Raku Son; raku.son@riken.jp;
### Reviewed by:
  Akiko Oguchi; akiko.oguchi@riken.jp;
  Sho Sekito; sekito.shou.25u@st.kyoto-u.ac.jp
  Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
  Yasuhiro Murakawa; yasuhiro.murakawa@riken.jp, murakawa.yasuhiro.0r@kyoto-u.ac.jp

### The following log file uses test files named as Test_S1_L001_R1_001.fastq.gz and Test_S1_L001_R2_001.fastq.gz

### Requirements
Cellranger v7.0.0 (https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
STAR v2.7.10a (https://github.com/alexdobin/STAR)
umi-tools v1.1.2 (https://github.com/CGATOxford/UMI-tools/releases/tag/1.1.2)
samtools v1.15.1 (http://www.htslib.org/)
bedtools v2.30.0 (https://github.com/arq5x/bedtools2)

### 0) Fastq files are in the working directory
ls
### Test_S1_L001_R1_001.fastq.gz  Test_S1_L001_R2_001.fastq.gz
\

### Move fastqs in each directory to process multiple samples at one time
for infile in *R1_001.fastq.gz; do base=$(basename ${infile} _R1_001.fastq.gz); mkdir ${base}_fastq; mv ${infile} ${base}_fastq/; mv ${base}_R2_001.fastq.gz  ${base}_fastq/; done

ls
### Test_S1_L001_fastq

### 1) Cellranger count
### Cell barcodes and count data from Cellranger are used in the subsequent analysis
### Although STARsolo also filters cell barcodes and produces count data, the result usually contains less #cells compared with Cellranger, which sometimes misses specific clusters

export PATH=~/cellranger-7.0.0:$PATH
for infile in *_fastq;
do base=$(basename ${infile} _fastq);
cellranger count --id=${base}_Cellranger \
--fastqs=${infile} \
--transcriptome=/local/home/ubuntu/refdata-gex-GRCh38-2020-A \
--localmem=256;
done

### Prepare whitelist of cell barcodes
mkdir whitelists
for infile in *_fastq;
do base=$(basename ${infile} _fastq);
zcat ${base}_Cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed -e 's/-1//g' > whitelists/${base}_whitelist.txt;
done

### 2) STAR mapping
### Without --limitBAMsortRAM option, STAR sometimes outputs error as
### EXITING because of fatal ERROR: not enough memory for BAM sorting:
### SOLUTION: re-run STAR with at least --limitBAMsortRAM (some specific number)
### 60000000000 is usually fine but in some deeply sequenced samples, 120000000000 was required

### Usually, allowing to open a large number of files is needed
### Otherwise, STAR will outputs error as
### BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file Test_S1_L001_STARsolo/Test_S1_L001__STARtmp//BAMsort/19/45
### SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
ulimit -n 10000

for infile in *_fastq;
do base=$(basename ${infile} _fastq);
STAR \
--runThreadN 32 \
--genomeDir /local/home/ubuntu/ref/reference_20220827/human_GRCh38_gencodev41/GRCh38_index/ \
--readFilesIn ${base}_fastq/${base}_R1_001.fastq.gz ${base}_fastq/${base}_R2_001.fastq.gz \
--soloCBwhitelist whitelists/${base}_whitelist.txt \
--soloBarcodeMate 1 --clip5pNbases 39 0 \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 16 \
--soloUMIstart 17 --soloUMIlen 10 \
--soloStrand Reverse \
--outFileNamePrefix ${base}_STARsolo/${base}_ \
--outSAMtype BAM SortedByCoordinate \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIdedup 1MM_Directional_UMItools \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--limitBAMsortRAM 60000000000;
done

### Oct 22 10:30:47 ..... started STAR run
### Oct 22 10:30:47 ..... loading genome
### Oct 22 10:43:11 ..... started mapping
### Oct 22 10:47:13 ..... finished mapping
### Oct 22 10:47:33 ..... started Solo counting
### Oct 22 10:47:43 ..... finished Solo counting
### Oct 22 10:47:43 ..... started sorting BAM
### Oct 22 10:49:11 ..... finished successfully

### 3) Extract reads that start with unencoded-G
### Extract uniquely mapped reads (Read1 only)
for infile in *_fastq;
do base=$(basename ${infile} _fastq);
samtools view -@ 32 -hbf 64 -q 255 ${base}_STARsolo/${base}_Aligned.sortedByCoord.out.bam > ${base}_STARsolo/${base}_Aligned.sortedByCoord.out_unique_R1.bam
done

mkdir Uniquely_mapped_R1
mv *_STARsolo/*_Aligned.sortedByCoord.out_unique_R1.bam Uniquely_mapped_R1/
cd Uniquely_mapped_R1/

### Extract unencoded-G reads
bash ../STARsolo_SoftclipG_221224.sh -t 32
### [bam_sort_core] merging from 0 files and 32 in-memory blocks...
cd ..

### 4) Deduplicate by umi-tools
### Prepare index files for umi-tools
for infile in *_fastq;
do base=$(basename ${infile} _fastq);
samtools index -@ 32 Uniquely_mapped_R1/SoftclipG_${base}_Aligned.sortedByCoord.out_unique_R1.bam;
done

### Umi-tools is used without --output-stats option since it sometimes produces error due to larger memory usage
### Since umi-tools only uses single-thread, this process is recommended to be parallelized as needed to save the total process time for multiple samples

for infile in *_fastq;
do base=$(basename ${infile} _fastq);
umi_tools dedup --per-cell \
-I Uniquely_mapped_R1/SoftclipG_${base}_Aligned.sortedByCoord.out_unique_R1.bam \
--extract-umi-method=tag --umi-tag=UR --cell-tag=CR \
-S Uniquely_mapped_R1/SoftclipG_${base}_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR.bam;
done

### 5) Filter by valid cell barcodes
### Modify whitelist to fit samtools input
mkdir cell_barcodes
for infile in *_fastq;
do base=$(basename ${infile} _fastq);
awk '{print "CB:Z:"$1}' whitelists/${base}_whitelist.txt > cell_barcodes/${base}_cell_barcode.txt;
done

for infile in *_fastq;
do base=$(basename ${infile} _fastq);
cd Uniquely_mapped_R1/
samtools view -@ 32 -H SoftclipG_${base}_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR.bam > SAM_header
samtools view -@ 32 SoftclipG_${base}_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR.bam | LC_ALL=C grep -F -f ../cell_barcodes/${base}_cell_barcode.txt > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -@ 32 -b filtered.sam > SoftclipG_${base}_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.bam
rm SAM_header
rm filtered_SAM_body
rm filtered.sam
cd ..;
done

mkdir SoftclipG_deduplicated_filtered_bam/
mv Uniquely_mapped_R1/SoftclipG_*_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.bam SoftclipG_deduplicated_filtered_bam/

### 6) Convert bam file to 1-base 5'-end tag file and performs QC
cd  SoftclipG_deduplicated_filtered_bam/

### Extract cell barcodes from softclipG bam file in the same order as in bam files
for infile in *bam;
 do
   base=$(basename ${infile} .bam)
   samtools view -@ 32 ${infile} \
 | awk 'BEGIN{OFS="\t"}{print $23}' \
 > ${base}_cell_barcode_tmp.txt
done

### Convert bam files to 5'-end bed files (CTSS.bed) with cell barcodes
### This part is separated in an individual script since bedtools only uses single-thread and it is recommended to be parallelized as needed to save the total process time for multiple samples
bash ../STARsolo_Cell_barcode_CTSS_bed_20221123.sh

### Remove tmp file
rm *_cell_barcode_tmp.txt

### We routinely count 5'-end tags on FANTOM5 promoter and FANTOM-NET enhancer regions as quality check
### Also 5'-end bed files can be visualized on IGV or UCSC by converting them to CPM-normalized bigwig
bash ../STARsolo_Counts_CTSS_bed_CPM_bigwig_240119.sh \
-G ~/ref/reference_20220827/human_GRCh38_gencodev41/GRCh38.primary_assembly.chrom.sizes \
-t 32 \
-p ~/ref/reference_20220703/human_GRCh38_gencodev40/hg38_fair+new_CAGE_peaks_phase1and2.6.bed \
-e ~/ref/reference_20220703/human_GRCh38_gencodev40/Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.6.bed

### 7) Perform enhancer call
### Enhancer call is based on Andersson's script (https://github.com/anderssonrobin/enhancers) with small modifications. (enhancers/scripts/fixed_bidir_enhancers_10bp)
### We use TSS±300bp of all transcripts except lncRNAs as masking files
### This filtering can be disabled by removing -m option in the script

### 7-1) Enhancer call from all cells
mkdir ../enhancer_call

### Prepare path of CTSS,bed files and perform enhancer call
for infile in *.CTSS.fwd.rev.cell.barcode.bed;
do base=$(basename ${infile} .CTSS.fwd.rev.cell.barcode.bed);
find `pwd` -name ${infile} > ../enhancer_call/bedlist_${base}.txt;
../enhancers/scripts/fixed_bidir_enhancers_10bp \
-f ../enhancer_call/bedlist_${base}.txt \
-m ~/ref/reference_20220827/human_GRCh38_gencodev41/gencodev41_masking_exceptlncRNA/Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.mask.6.bed \
-o ../enhancer_call/enhancercall_${base}/;
done

cd ../enhancer_call/enhancercall_SoftclipG_Test_S1_L001_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered
cd ..
cd ..

### 7-2) Enhancer call from specific cell clusters
### The same can be performed on specific cell clusters by subsetting 5'-end bed files
### List of cell barcodes needs to be specified by lists
### Here we are showing the example using cell barcode list of Cluster0_cell_barcode.txt
mkdir enhancer_call_cluster0
cd enhancer_call_cluster0

grep ../SoftclipG_deduplicated_filtered_bam/SoftclipG_Test_S1_L001_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.CTSS.fwd.rev.cell.barcode.bed -f ../Cluster0_cell_barcode.txt > SoftclipG_Test_S1_L001_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.CTSS.fwd.rev.cell.barcode_cluster0.bed

find `pwd` -name SoftclipG_Test_S1_L001_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.CTSS.fwd.rev.cell.barcode_cluster0.bed > bedlist_cluster0.txt

../enhancers/scripts/fixed_bidir_enhancers_10bp -f bedlist_cluster0.txt -m ~/ref/reference_20220827/human_GRCh38_gencodev41/gencodev41_masking_exceptlncRNA/Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.mask.6.bed -o enhancer_call_cluster0

cd enhancer_call_cluster0/

### Compress bed files when the analysis is completed
cd ..
pigz -p 32 SoftclipG_deduplicated_filtered_bam/*.CTSS.fwd.rev.cell.barcode.bed
