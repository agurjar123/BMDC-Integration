# Data Directory

Place your Seurat-exported count matrices here.

## Required Files

For each sample (numbered 1, 2, 3):
- `counts_N.mtx` - Count matrix in Matrix Market format
- `features_N.tsv` - Gene features (gene IDs and symbols)
- `barcodes_N.tsv` - Cell barcodes
- `metadata_N.csv` - Cell metadata

## Example Structure

```
counts/
├── counts_1.mtx
├── features_1.tsv
├── barcodes_1.tsv
├── metadata_1.csv
├── counts_2.mtx
├── features_2.tsv
├── barcodes_2.tsv
├── metadata_2.csv
├── counts_3.mtx
├── features_3.tsv
├── barcodes_3.tsv
└── metadata_3.csv
```

## Exporting from Seurat

Use the following R code to export your Seurat object:

```r
library(Seurat)
library(Matrix)

# Export count matrix
writeMM(GetAssayData(seurat_obj, slot='counts'), 'counts_1.mtx')

# Export features
write.table(
  data.frame(rownames(seurat_obj)),
  'features_1.tsv',
  sep='\t',
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE
)

# Export barcodes
write.table(
  data.frame(colnames(seurat_obj)),
  'barcodes_1.tsv',
  sep='\t',
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE
)

# Export metadata
write.csv(seurat_obj@meta.data, 'metadata_1.csv')
```

