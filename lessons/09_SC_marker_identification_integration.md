---
title: "Single-cell RNA-seq: Marker identification"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Thursday, June 6,2019
---

Approximate time: 45 minutes

## Learning Objectives:

* Understand how to determine markers of individual clusters
* Understand the iterative processes of clustering and marker identification

# Single-cell RNA-seq marker identification

Now that we have identified our desired clusters, we can move on to marker identification, which will allow us to verify the identity of certain clusters and help surmise the identity of any unknown clusters. 

<img src="../img/sc_workflow_integration.png" width="800">

***

_**Goals:**_ 
 
 - _To **determine the gene markers** for each of the clusters_
 - _To **identify cell types** of each cluster using markers_
 - _To determine whether need to **re-cluster based on cell type markers**, perhaps clusters need to be merged or split_

_**Challenges:**_
 
 - _Over-interpretation of the results_
 - _Combining different types of marker identification_

_**Recommendations:**_
 
 - _Think of the results as hypotheses that need verification. Inflated p-values can lead to over-interpretation of results (essentially each cell is used as a replicate). Top markers are most trustworthy._
 _Identify all markers conserved between conditions for each cluster_
 - _Identify markers that are differentially expressed between specific clusters_

***

Remember that we had the following questions from the clustering analysis:

1. *What is the cell type identity of clusters 17 and 20?*
3. *Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?*
4. *Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?*

Remember that we have a few different types of marker identification, each with their own benefits and drawbacks:

1. **Identification of all markers for each cluster**
	
	- *Useful for identifying unknown clusters and improving confidence in hypothesized cell types.*

2. **Identification of conserved markers for each cluster regardless of condition** 
	
	- *Useful for identifying unknown clusters and improving confidence in cell type identities when more than one condition. Returns fewer, but higher confidence markers of cell type. Tends to take much longer to run.*  	

3. **Marker identification between specific clusters**
	
	- *Useful for determining differences in gene expression between clusters with markers that are similar in the above analyses.*

Since we have more than a single condition, we can explore all methods of marker identification.

## Identification of all markers for each cluster

This analysis compares each cluster against all others and outputs the genes that are differentially expressed/present using the `FindAllMarkers()` function. 

We have already accomplished this with the `ctrl` and `stim` samples separately. Let's use our knownledge to complete the exercises below.

***

**[Exercises](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html)**

1. **Find all markers** for the `combined` data using the `FindAllMarkers()` function with arguments to return only the positive markers and those markers with log2 fold-change greater than 0.25.
2. **Annotate** the markers with gene descriptions.
3. **Reorder the columns** to be in the order shown below.

	<p align="center">
	<img src="../img/all_comb_markers.png" width="800">
	</p>

4. **Arrange rows** by cluster, then by p-adjusted values
5. **Save** our rearranged marker analysis results to a file called `combined_all_markers.csv` in the `results` folder.
6. **Extract the top 5 markers** by log2 fold change for each cluster and view them.
7. Based on these marker results, **determine whether the markers make sense** for our hypothesized identities for each cluster:

	| Cell Type | Clusters |
	|:---:|:---:|
	| CD14+ Monocytes | 0 | 
	| FCGR3A+ Monocytes | 7 |
	| Conventional dendritic cells | 13 |
	| Plasmacytoid dendritic cells | 18 |
	| Naive B cells | 4, 15 |
	| Activated B cells | 14 |
	| Naive CD4+ T cells | 1, 2, 9, 16 |
	| Activated CD4+ T cells | 3 |
	| Naive CD8+ T cells| 11 |
	| Activated (cytotoxic) CD8+ T cells| 6, 10 |
	| NK cells | 5 |
	| Megakaryocytes | 12 |
	| Erythrocytes | 19 |
	| Stressed / dying cells | 8 |
	| Unknown | 17, 20 |

***

If there were any questions about the identity of any clusters, exploring the cluster's markers would be the first step. Let's look at the `ann_comb_markers`, filtering for cluster 17. ATAD2 can be associated with activation, so maybe cluster 17 corresponds to activated CD4+ T cells.

<p align="center">
<img src="../img/sc_integ_cluster17.png" width="800">
</p>

We also had questions regarding the identity of cluster 20. We can look at the markers of cluster 20 to try to resolve the identity:

<p align="center">
<img src="../img/sc_integ_cluster20.png" width="800">
</p>

The enriched genes appear to be markers of proliferation. In addition, we see expression of genes encoding T cell inhibitory receptors. However, cluster 20 had low levels of the T cell markers CD3D, CD4 and CD8. So this cluster is still a bit of a mystery.


## Identification of conserved markers in all conditions

Identifying conserved markers allows for identifying only those genes the are significantly differentially expressed relative to the other clusters for all conditions. This function is most useful to run if unsure of the identity for a cluster after running the `FindAllMarkers()`. You could run it on all clusters if you wanted to, but it takes a while to run, so we are just going to run it on the unknown clusters 17 and 20.

The function we will use is the `FindConservedMarkers()`, which has the following structure:

**`FindConservedMarkers()` syntax:**

```r
FindConservedMarkers(seurat_obj,
                     ident.1 = cluster,
                     grouping.var = "group",
                     only.pos = TRUE)
```

The function **accepts a single cluster at a time**, so if we want to have the function run on all clusters, then we can use the `map` family of functions to iterate across clusters. 

Since these functions will **remove our row names** (gene names), we need to transfer the row names to columns before mapping across clusters. We also need a column specifying **to which cluster the significant genes correspond**.

To do that we will create our own function to:

1. Run the `FindConservedMarkers()` function
2. Transfer row names to a column using `rownames_to_column()` function
3. Create the column of cluster IDs using the `cbind()` function

```r
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
        FindConservedMarkers(combined,
                             ident.1 = cluster,
                             grouping.var = "sample",
                             only.pos = TRUE) %>%
                rownames_to_column(var = "gene") %>%
                cbind(cluster_id = cluster, .)
}
```

Since we want the output of the `map` family of functions to be a **dataframe with each cluster output bound together by rows**, we will use the `map_dfr()` function.

Remember the map family of functions uses the following syntax:

**`map` family syntax:**

```r
map_dfr(inputs_to_function, name_of_function)
```

Now, let's find the conserved markers for clusters 17 and 20. 

```r
# Iterate function across desired clusters
conserved_markers <- map_dfr(c(17,20), get_conserved)
```

> **NOTE:** If you wanted to run this on all clusters, you could input `0:20` instead of `c(17,20)`; however, it would take quite a while to run.

To better analyze the output, we can include the gene descriptions as well.

```r
# Extract the gene descriptions for each gene
gene_descriptions <- unique(annotations[, c("gene_name", "description")])

# Merge gene descriptions with markers
ann_conserved_markers <- left_join(x = conserved_markers,
                                   y = gene_descriptions,
                                   by = c("gene" = "gene_name"))
```

<p align="center">
<img src="../img/sc_integ_marker_unknown.png" width="800">
</p>

For clusters 17 and 20, we see many of the conserved enriched genes encode inhibitory receptors, such as TIGIT and LAG3, which can be indicative of exhausted T cells.

## Identifying gene markers for each cluster

The last set of questions we had regarding the analysis involved whether the clusters corresponding to the same cell types have biologically meaningful differences. Sometimes the list of markers returned don't sufficiently separate some of the clusters. 

Again, we performed this analysis with the single samples, so we are going to perform it by completing the exercises below.

***

**Exercises**

1. Determine differentiating markers for CD8+ T cells - clusters 6 versus 10 - using the `FindMarkers()` function.
2. **Annotate** the markers with gene descriptions.
3. **Reorder the columns** to be in the order shown below.

	<p align="center">
	<img src="../img/sc_cd8_t_markers_diff.png" width="800">
	</p>

4. **Arrange rows** by `avg_logFC` values
5. **Save** our rearranged marker analysis results to a file called `cluster6vs10_markers.csv` in the `results` folder.
6. Based on these marker results, **determine whether we need to separate** clusters 6 and 10 as their own clusters.
7. **Extra credit:** Repeat above steps for the clusters assigned to `Naive CD4+ T cells`, in addition to repeating for `Naive B cells`.

***


Now taking in all of this information, we can surmise the cell types of the different clusters and plot the cells with cell type labels. 

| Cell Type | Clusters |
|:---:|:---:|
| CD14+ Monocytes | 0 | 
| FCGR3A+ Monocytes | 7 |
| Conventional dendritic cells | 13 |
| Plasmacytoid dendritic cells | 18 |
| Naive B cells | 4, 15 |
| Activated B cells | 14 |
| Naive CD4+ T cells | 1, 2, 9, 16 |
| Activated CD4+ T cells | 3 |
| Naive CD8+ T cells| 11 |
| Activated (cytotoxic) CD8+ T cells| 6, 10 |
| NK cells | 5 |
| Megakaryocytes | 12 |
| Erythrocytes | 19 |
| Stressed/dying cells | 8 |
| Exhausted T cells | 17
| Proliferating unknown | 20 |

We can then reassign the identity of the clusters to these cell types:

```r
combined <- RenameIdents(object = combined,
                         "0" = "CD14+ monocytes",
                         "1" = "Naive CD4+ T cells",
                         "2" = "Naive CD4+ T cells",
                         "3" = "Activated CD4+ T cells",
                         "4" = "Naive B cells",
                         "5" = "NK cells",
                         "6" = "Activated (cytotoxic) CD8+ T cells",
                         "7" = "FCGR3A+ Monocytes",
                         "8" = "Stressed/dying cells",
                         "9" = "Naive CD4+ T cells",
                         "10" = "Activated (cytotoxic) CD8+ T cells",
                         "11" = "Naive CD8+ T cells",
                         "12" = "Megakaryocytes",
                         "13" = "Conventional dendritic cells",
                         "14" = "Activated B cells",
                         "15" = "Naive B cells",
                         "16" = "Naive CD4+ T cells",
                         "17" = "Exhausted T cells",
                         "18" = "Plasmacytoid dendritic cells",
                         "19" = "Erythrocytes",
                         "20" = "Proliferating unknown")

DimPlot(object = combined, 
        reduction = "umap", 
        label = TRUE,
        pt.size = 0.5,
        repel = T,
        label.size = 6)
```

<p align="center">
<img src="../img/sc_integ_umap_labelled.png" width="800">
</p>

***

**Exercises**

1. Remove the stressed cells using the `subset()` function.
2. Visualize the clusters using `DimPlot()`.
3. Use the `write_rds()` function to save the final labelled `combined` object to the `results` folder, called `combined_labelled_res0.8.rds`.

***

Now that we have our clusters defined and the markers for each of our clusters, we have a few different options:

- Experimentally validate intriguing markers for our identified cell types.
- Perform differential expression analysis between conditions `ctrl` and `stim`
	- Biological replicates are **necessary** to proceed with this analysis.
- Trajectory analysis, or lineage tracing, could be performed if trying to determine the progression between cell types or cell states. For example, we could explore any of the following using this type of analysis:
	- Differentiation processes
	- Expression changes over time
	- Cell state changes in expression

***

[Explore the answer key for the exercises in this lesson]()

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
