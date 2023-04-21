## clean spatial ATAC-seq data
Latch workflow for remediating lane artifacts in spatial ATAC-seq experiments;
specifically for data generated via DBiT-seq (see [Deng, 2022](https://www.nature.com/articles/s41586-022-05094-1)).

Workflow takes the following outputs from Cellragnger ATAC,
- fragments.tsv.gz
- singlecell.csv

and returns a 'cleaned' fragments.tsv.gz where hot row/columns are downsampled to be 
within x standard deviations of the mean of row/column medians.  

Output table is sorted for continuous chromosome blocks, and compressed with gbzip for ArchR.

Output can be used for analysis with ArchR, Seurat, and other scATAC-seq packages.

