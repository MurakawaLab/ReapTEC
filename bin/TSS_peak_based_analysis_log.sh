cd ~/TC_analysis/5scCTSSbed_All

# *_Aligned.sortedByCoord.out_unique_deduplicated_UR_CR_filtered.CTSS.fwd.rev.cell.barcode.bed
# files should be first sorted by "LANG=C sort -k 4"
# Please sort ctss.bed files before running this script.

for infile in *_sorted.CTSS.bed; do
echo ${infile}
base=$(basename ${infile} _sorted.CTSS.bed)
bedtools intersect -wa -wb -s \
-a ~/TC_analysis/20230928_TC_anno_Gene_rank_bed6.bed \
-b ${infile} \
| cut -f4,10,11 > ${base}.txt

cat ${base}.txt | cut -f1 | sort | uniq | cat -n | LANG=C sort -k 2 > ${base}.feature.number.sort.txt
cat ${base}.feature.number.sort.txt | cut -f2 > ${base}.features.tsv
cat ${base}.txt | cut -f2 | sort | uniq | cat -n | LANG=C sort -k 2 > ${base}.barcode.number.sort.txt
cat ${base}.barcode.number.sort.txt | cut -f2 > ${base}.barcodes.tsv

LANG=C sort -k 1 ${base}.txt > ${base}.feature.sort.txt
join -o '1.1,0,2.2,2.3' -1 2 -2 1 ${base}.feature.number.sort.txt ${base}.feature.sort.txt > ${base}.featureNo.txt
LANG=C sort -k 3 ${base}.featureNo.txt > ${base}.featureNo.barcode.sort.txt
join -o '1.1,2.1,1.2,0,1.4' -1 3 -2 2 ${base}.featureNo.barcode.sort.txt ${base}.barcode.number.sort.txt > ${base}.featureNo.barcodeNo.txt
cat ${base}.featureNo.barcodeNo.txt | cut -f1,2,5 -d' ' > ${base}.featureNo.barcodeNo.trim.txt
cat ${base}.featureNo.barcodeNo.trim.txt \
| awk 'BEGIN {OFS="\t"} {if ($3 != 0) {print $0}}' \
> ${base}.featureNo.barcodeNo.trim_trim0.txt

###
feature_no=$(wc -l ${base}.features.tsv | cut -f1 -d " ")
barcode_no=$(wc -l ${base}.barcodes.tsv | cut -f1 -d " ")
read_count=$(wc -l ${base}.featureNo.barcodeNo.trim_trim0.txt | cut -f1 -d " ")

echo %%MatrixMarket matrix coordinate real general >> ${base}.tmp.txt
echo ${feature_no}$'\t'${barcode_no}$'\t'${read_count} >> ${base}.tmp.txt
cat ${base}.tmp.txt ${base}.featureNo.barcodeNo.trim_trim0.txt > ${base}.matrix.mtx

rm ${base}.txt
rm ${base}.feature.number.sort.txt
rm ${base}.barcode.number.sort.txt
rm ${base}.feature.sort.txt
rm ${base}.featureNo.txt
rm ${base}.featureNo.barcodeNo.txt
rm ${base}.featureNo.barcode.sort.txt
rm ${base}.tmp.txt
rm ${base}.featureNo.barcodeNo.trim.txt
rm ${base}.featureNo.barcodeNo.trim_trim0.txt

mkdir ${base}

mv ${base}.features.tsv features.tsv
sed "s/CB:Z://g" ${base}.barcodes.tsv > barcodes.tsv
rm ${base}.barcodes.tsv
mv ${base}.matrix.mtx matrix.mtx
gzip *tsv
gzip *mtx
mv *gz ./${base}
done

#An example of R code to show how to inport the count data created above to Seurat object. 
matrix_dir = "~/TC_analysis/5scCTSSbed_All/20240411_Test/SoftclipG_D1_Bulk_5sc_1/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
D1_Bulk_5sc_1.data <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(D1_Bulk_5sc_1.data) = barcode.names$V1
rownames(D1_Bulk_5sc_1.data) = feature.names$V1
D1_Bulk_5sc_1 = CreateSeuratObject(counts = Matrix::Matrix(as.matrix(D1_Bulk_5sc_1.data),sparse = T), project = "D1_Bulk_5sc_1")

# CD4bulk_5sc_TCver_TCRremoved.removeMono.harmonyintegration_FigS14.rds 
# is a file containing Seurat object used in Supplementary figure S14 in the manuscript (Oguchi et al. Science).
# The data from 5′ scRNA-seq of bulk CD4+ T cells were integrated using Harmony.