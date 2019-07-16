### Creating count data object

Generally, all single-cell RNA-seq datasets, regardless of technology or pipeline, will contain **three files**:

1. a file with the **gene IDs**, representing all genes quantified
2. a file with the **cell IDs**, representing all cells quantified
3. a **matrix of counts** per gene for every cell


We can explore these files by clicking on the `data/ctrl_raw_feature_bc_matrix` folder:

- **`barcodes.tsv`:** cellular barcodes present in dataset

  <p align="center">
  <img src="../img/cell_ids_new.png" width="180">
  </p>
  
- **`genes.tsv`:** IDs of quantified genes

  <p align="center">
  <img src="../img/genes.png" width="300">
  </p>

- **`matrix.mtx`:** a matrix of count values, where rows are associated with the gene IDs above and columns correspond to the cellular barcodes. Note that there are many zero values in this matrix.

  <p align="center">
  <img src="../img/sparse_matrix.png">
  </p>

We can create a count matrix using these files. However, instead of creating a standard count matrix, we will create a **sparse matrix** to improve the amount of space, memory and CPU required to work with our huge count matrix. 

We will use `readMM()` function from the **Matrix** package to turn our standard matrix into a sparse matrix. The `genes.tsv` file should correspond to the genes or row names of the matrix, while `barcodes.tsv` corresponds to the cells or columns.

```r
# Read in `matrix.mtx`
counts <- readMM("data/ctrl_raw_feature_bc_matrix/matrix.mtx")

# Read in `genes.tsv`
genes <- read_tsv("data/ctrl_raw_feature_bc_matrix/genes.tsv", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cell_ids <- read_tsv("data/ctrl_raw_feature_bc_matrix/barcodes.tsv", col_names = FALSE)$X1
```

Then we can add row names to the count matrix to be the gene IDs and the column names of the count matrix to be the cell IDs.

```r
# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
```

We could use this data for downstream QC analysis. However, this would take a long time if we had multiple samples. A quicker way to load multiple samples is to use the Seurat R package, which has a specific function for reading in 10X data, called `read10X()`. 

> **NOTE:** If using other droplet-based methods for library preparation, the above method would be needed to perform the QC. We have [additional materials](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html) available based on creation of the count matrix in this way.

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
