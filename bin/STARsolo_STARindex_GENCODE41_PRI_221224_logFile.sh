#!/usr/bin/env bash

### GENCODE v41
### Primary assembly (PRI) 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz

### Comprehensive PRI gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz

gunzip *

### Make GRCh38 PRI GENCODE v41 index
STAR --runThreadN 20 \
 --runMode genomeGenerate \
 --genomeDir ~/ref/reference_20220827/human_GRCh38_gencodev41/GRCh38_index \
 --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
 --sjdbGTFfile gencode.v41.primary_assembly.annotation.gtf \
 --sjdbOverhang 149

# Aug 27 02:03:42 ..... started STAR run
# Aug 27 02:03:42 ... starting to generate Genome files
# Aug 27 02:04:27 ..... processing annotations GTF
# Aug 27 02:04:59 ... starting to sort Suffix Array. This may take a long time...
# Aug 27 02:05:14 ... sorting Suffix Array chunks and saving them to disk...
# Aug 27 02:25:17 ... loading chunks from disk, packing SA...
# Aug 27 02:26:46 ... finished generating suffix array
# Aug 27 02:26:46 ... generating Suffix Array index
# Aug 27 02:33:19 ... completed Suffix Array index
# Aug 27 02:33:20 ..... inserting junctions into the genome indices
# Aug 27 02:38:37 ... writing Genome to disk ...
# Aug 27 02:38:39 ... writing Suffix Array to disk ...
# Aug 27 02:39:13 ... writing SAindex to disk
# Aug 27 02:39:17 ..... finished successfully