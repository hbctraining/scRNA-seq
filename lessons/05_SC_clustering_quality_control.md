---
title: "Single-cell RNA-seq: Clustering Analysis"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Thursday, June 6, 2019
---

Approximate time: 90 minutes

## Learning Objectives:

* Evaluate whether clustering artifacts are present 
* Determine the quality of clustering with PCA, tSNE and UMAP plots and understand when to re-cluster
* Assess known cell type markers to hypothesize cell type identities of clusters

# Single-cell RNA-seq clustering analysis


Now that we have our high quality cells, we want to know the different cell types present within our population of cells. 

<img src="../img/sc_workflow.png" width="800">

***

_**Goals:**_ 
 
 - _To **generate cell type-specific clusters** and use known markers to determine the identities of the clusters._
 - _To **determine whether clusters represent true cell types or cluster due to biological or technical variation**, such as clusters of cells in the S phase of the cell cycle, clusters of specific batches, or cells with high mitochondrial content._

_**Challenges:**_
 
 - _Clustering so that **cells of the same cell type from different conditions cluster together**_
 - _**Removing unwanted variation** so that we do not have cells clustering by artifacts_
 - _**Identifying the cell types** of each cluster_
 - _Maintaining patience as this can be a highly iterative process between clustering and marker identification (sometimes even going back to the QC filtering)_

_**Recommendations:**_
 
 - _Have a good idea of your expectations for the **cell types to be present** prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating_
 - _If you have **more than one condition**, it's often helpful to perform integration to align the cells_
 - _**Regress out** number of UMIs, mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering_
 - _Identify any junk clusters for removal. Possible junk clusters could include those with high **mitochondrial content** and low UMIs/genes_
 - _If **not detecting all cell types as separate clusters**, try changing the resolution or the number of PCs used for clustering_

***

## Exploration of quality control metrics

To determine whether our clusters might be due to artifacts such as cell cycle phase or mitochondrial expression, it can be useful to explore these metrics visually to see if any clusters exhibit enrichment or are different from the other clusters. However, if enrichment or differences are observed for particular clusters it may not be worrisome if it can be explained by the cell type. 

We can start by exploring the distribution of cells per cluster:

```{r cell_counts}
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_control, 
                     vars = c("ident")) %>% 
        dplyr::count(ident) %>% 
        spread(ident, n)

# View table
View(n_cells)
```

<p align="center">
<img src="../img/SC_clustercells_loadObj.png" width="800">
</p>

Tho acquire the different cluster QC metrics, we can use the `FetchData()` function from Seurat, perform some data wrangling, and plot the metrics with ggplot2. We will start by exploring the distribution of cells in each sample and in the different phases of the cell cycle to view by UMAP and PCA.

First we will acquire the cell cycle and UMAP coordinate information to view by UMAP:

```{r plot_feature_tsne, fig.width=10, fig.height=5}
# Establishing groups to color plots by
group_by <- c("Phase")

# Getting coordinates for cells to use for UMAP and associated grouping variable information
class_umap_data <- FetchData(seurat_control, 
                             vars = c("ident", "UMAP_1", "UMAP_2", group_by))

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_control, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
        as.data.frame() %>% 
        group_by(ident) %>%
        summarise(x=mean(UMAP_1), y=mean(UMAP_2))
```

> **NOTE:** How did we know in the `FetchData()` function to include `UMAP_1` to obtain the UMAP coordinates? The [Seurat cheatsheet](https://satijalab.org/seurat/essential_commands.html) describes the function as being able to pull any data from the expression matrices, cell embeddings, or metadata. 
> 
> For instance, if you explore the `seurat_control@reductions` list object, the first component is for PCA, and includes a slot for `cell.embeddings`. We can use the column names (`PC_1`, `PC_2`, `PC_3`, etc.) to pull out the coordinates or PC scores corresponding to each cell for each of the PCs. 
> 
> We could do the same thing for UMAP:
> 
> ```r
> # Extract the UMAP coordinates for the first 10 cells
> seurat_control@reductions$umap@cell.embeddings[1:10, 1:2]
>```
>
> The `FetchData()` function just allows us to extract the data more easily.


In addition, we can acquire the same metrics to view by PCA:

```r
# Getting coordinates for cells to use for PCA and associated grouping variable information
class_pca_data <- FetchData(seurat_control, 
                            vars = c("ident", "PC_1", "PC_2", group_by))

# Adding cluster label to center of cluster on PCA
pca_label <- FetchData(seurat_control, 
                       vars = c("ident", "PC_1", "PC_2"))  %>%
        as.data.frame() %>%
        mutate(ident = seurat_control@active.ident) %>%
        group_by(ident) %>%
        summarise(x=mean(PC_1), y=mean(PC_2))
```

Then, we can plot the cell cycle by UMAP and PCA:

```r
# Function to plot UMAP and PCA as grids
  map(group_by, function(metric) {
    plot_grid(
      ggplot(class_umap_data, aes(UMAP_1, UMAP_2)) +
        geom_point(aes_string(color = metric), alpha = 0.7) +
        scale_color_brewer(palette = "Set2")  +
        geom_text(data=umap_label, aes(label=ident, x, y)),
      ggplot(class_pca_data, aes(PC_1, PC_2)) +
        geom_point(aes_string(color = metric), alpha = 0.7) +
        scale_color_brewer(palette = "Set2")  +
        geom_text(data=pca_label, 
                  aes(label=ident, x, y)),
      nrow = 1, 
      align = "v")
  })

```

<p align="center">
<img src="../img/SC_phase_umap_pca_loadObj.png" width="800">
</p>

Next we will explore additional metrics, such as the number of UMIs and genes per cell, S-phase and G2M-phase markers, and mitochondrial gene expression by UMAP. 

```{r dim_features, fig.width=10, fig.height=5}
# Determine metrics to plot present in seurat_control@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

# Extract the UMAP coordinates for each cell and include information about the metrics to plot
qc_data <- FetchData(seurat_control, 
                     vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))

# Plot a UMAP plot for each metric
map(metrics, function(qc){
        ggplot(qc_data,
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=qc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(qc)
}) %>%
        plot_grid(plotlist = .)
```

<p align="center">
<img src="../img/SC_metrics_umpa_loadObj.png" width="800">
</p>

The metrics seem to be relatively even across the clusters, with the exception of the `nUMIs` and `nGene` exhibiting higher values in clusters 8 and 12 (and 11 to some extent). We will keep an eye on these clusters to see whether the cell types may explain the increase.

We can also explore how well our clusters separate by the different PCs; we hope that the defined PCs separate the cell types well. In the UMAP plots below, the cells are colored by their PC score for each respective principal component.

```{r feature_pcs, fig.width=10, fig.height=10}
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:14),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_control, 
                     vars = columns)

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:14), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)
```

<p align="center">
<img src="../img/SC_clusters_PCs_loadObj.png" width="800">
</p>

We can see how the clusters are represented by the different PCs. For instance, the genes driving `PC_2` exhibit higher expression in clusters 4 and 13, while clusters 6,7, and 14 show lower expression. We could look back at our genes driving this PC to get an idea of what the cell types might be:

```r
# Examine PCA results 
print(seurat_control[["pca"]], dims = 1:5, nfeatures = 5)

```

<p align="center">
<img src="../img/PC_print_loadObj.png" width="400">
</p>

With the CD79A gene and the HLA genes as positive markers of `PC_2`, we can hypothesize that clusters 4 and 13 correspond to B cells. This just hints at what the clusters identity could be, with the identities of the clusters being determined through a combination of the PCs. 

To truly determine the identity of the clusters and whether the `resolution` is appropriate, it is helpful to explore a handful of known markers for the cell types expected. 

## Exploring known cell type markers

With the cells clustered, we can explore the cell type identities by looking for known markers. The UMAP plot with clusters marked is shown, followed by the different cell types expected.

```{r UMAP_ref}
DimPlot(object = seurat_control, 
        reduction = "umap", 
        label = TRUE)
```

<p align="center">
<img src="../img/SC_umap_loadimg.png" width="800">
</p>


The `FeaturePlot()` function from seurat makes it easy to visualize a handful of genes using the gene IDs stored in the Seurat object. For example if we were interested in exploring known immune cell markers, such as:

| Cell Type | Marker |
|:---:|:---:|
| CD14+ monocytes | CD14, LYZ | 
| FCGR3A+ monocytes | FCGR3A, MS4A7 |
| Conventional dendritic cells | FCER1A, CST3 |
| Plasmacytoid dendritic cells | IL3RA, GZMB, SERPINF1, ITM2C |
| B cells | CD79A, MS4A1 |
| T cells | CD3D |
| CD4+ T cells | CD3D, IL7R, CCR7 |
| CD8+ T cells| CD3D, CD8A |
| NK cells | GNLY, NKG7 |
| Megakaryocytes | PPBP |
| Erythrocytes | HBB, HBA2 |

Seurat's `FeaturePlot()` function let's us easily explore the known markers on top of our t-SNE or UMAP visualizations. Let's go through and determine the identities of the clusters.

**CD14+ monocyte markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CD14", "LYZ"))
```

<p align="center">
<img src="../img/markers_CD14_monocytes.png" width="800">
</p>

CD14+ monocytes appear to correspond to clusters 0 and 5, and to a lesser extent, cluster 11.

**FCGR3A+ monocyte markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"))
```

<p align="center">
<img src="../img/markers_FCGR3A_monocytes.png" width="800">
</p>

FCGR3A+ monocytes markers distinctly highlight cluster 11.

**Conventional dendritic cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"))
```

<p align="center">
<img src="../img/markers_DCs.png" width="800">
</p>

The markers corresponding to conventional dendritic cells identify most of cluster 10.

**Plasmacytoid dendritic cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"))
```

<p align="center">
<img src="../img/markers_pDCs.png" width="800">
</p>

Plasmacytoid dendritic cells (pDCs) correspond to the part of cluster 10 that didn't mark the conventional dendritic cells (cDCs). This indicates that we may need to increase our `resolution` clustering parameter to separate out our pDCs from our cDCs. 

We could test out different resolutions by running the following code:

```r
# Assign identity of clusters
Idents(object = seurat_control) <- "RNA_snn_res.1.2"

# Plot the UMAP
DimPlot(seurat_control,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")

FeaturePlot(seurat_control, 
        reduction = "umap", 
        features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"))
```

Then, we could work our way back through the analysis, starting with the **`Exploration of quality control metrics`** section. However, in the interest of time, we will continue with the clusters we saw earlier, and go back to the resolution of "0.8".

```r
# Assign identity of clusters
Idents(object = seurat_control) <- "RNA_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_control,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")
```

**B cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CD79A", "MS4A1"))
```

<p align="center">
<img src="../img/markers_Bcells.png" width="800">
</p>

Clusters 4 and 13 have good expression of the B cell markers.

**T cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CD3D"))
```

<p align="center">
<img src="../img/markers_Tcells.png" width="600">
</p>

All T cells markers concentrate in clusters 1, 2, 3, 7, 8, 14, and 15.

**CD4+ T cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CD3D", "IL7R", "CCR7"))
```

<p align="center">
<img src="../img/markers_CD4Tcells.png" width="800">
</p>

The subset of T cells corresponding to the CD4+ T cells are clusters 1, 2, 3, 14, and 15.

**CD8+ T cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CD3D", "CD8A"))
```

<p align="center">
<img src="../img/markers_CD8Tcells.png" width="800">
</p>

Clusters 7 and 8 have high expression of the CD8+ T cell markers.

**NK cell markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("GNLY", "NKG7"))
```

<p align="center">
<img src="../img/markers_NKcells.png" width="800">
</p>

The NK cell markers are expressed in cluster 8, and to a lesser degree, cluster 7.

**Megakaryocyte markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("PPBP"))
```

<p align="center">
<img src="../img/markers_megakaryocytes.png" width="600">
</p>

The megakaryocyte markers seem to be expressed mainly in cluster 12.

**Erythrocyte markers**

```r
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("HBB", "HBA2"))
```

<p align="center">
<img src="../img/markers_erythrocytes.png" width="600">
</p>

There does not seem to be a cluster of erythrocytes, as the markers are spread around the different cell types. This is not a bad thing as the blood cells are often cell types to be excluded from the analysis, since they are not often informative about the condition of interest.

Based on these results, we can associate clusters with the cell types. However, we would like to perform a deeper analysis using marker identification before performing a final assignment of the clusters to a cell type.


| Cell Type | Clusters |
|:---:|:---:|
| CD14+ monocytes | 0, 5 | 
| FCGR3A+ monocytes | 11 |
| Conventional dendritic cells | 10 |
| Plasmacytoid dendritic cells | 10 |
| B cells | 4, 13 |
| T cells | 1, 2, 3, 7, 8, 14, 15 |
| CD4+ T cells | 1, 2, 3, 14, 15 |
| CD8+ T cells| 7, 8 |
| NK cells | 6, 7 |
| Megakaryocytes | 12 |
| Erythrocytes | - |
| Unknown | 9 |

> **NOTE:** With the pDCs, in addition to any other cluster that appears to contain two separate cell types, it's helpful to increase our clustering resolution to properly subset the clusters, as discussed above. Alternatively, if we still can't separate out the clusters using increased resolution, then it's possible that we had used too few principal components such that we are just not separating out these cell types of interest. To inform our choice of PCs, we could look at our PC gene expression overlapping the UMAP plots and determine whether our cell populations are separating by the PCs included.

Now we have a decent idea as to the cell types corresponding to the majority of the clusters, but some questions remain:

1. *What is the cell type identity of cluster 9?*
2. *Is cluster 7 a CD8+ T cell or an NK cell?*
3. *Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?*
4. *Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?*

Marker identification analysis can help us address all of these questions. The next step will be to perform marker identification analysis, which will output the genes that significantly differ in expression between clusters. Using these genes we can determine or improve confidence in the identities of the clusters/subclusters.

[Click here for next lesson]()

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
