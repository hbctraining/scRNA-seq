---
title: "Single-cell RNA-seq: Marker identification"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Monday, July 15, 2019
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
 - _To determine whether there's a need to **re-cluster based on cell type markers**, perhaps clusters need to be merged or split_

_**Challenges:**_
 
 - _Over-interpretation of the results_
 - _Combining different types of marker identification_

_**Recommendations:**_
 
 - _Think of the results as hypotheses that need verification. Inflated p-values can lead to over-interpretation of results (essentially each cell is used as a replicate). Top markers are most trustworthy._
 _Identify all markers conserved between conditions for each cluster_
 - _Identify markers that are differentially expressed between specific clusters_

***

Our clustering analysis resulted in the following clusters:

<p align="center">
<img src="../img/" width="800">
</p>

Remember that we had the following questions from the clustering analysis: (**These questions will change!**)

1. *What is the cell type identity of cluster 7?*
2. *Is cluster 6 a CD8+ T cell or an NK cell?* *Is cluster 13 a T cell or an NK cell?*
3. *Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?*
4. *Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?*

There are a few different types of marker identification that we can explore using Seurat to get to the answer of these questions. Each with their own benefits and drawbacks:

1. **Identification of all markers for each cluster:** this analysis compares each cluster against all others and outputs the genes that are differentially expressed/present. 
	- *Useful for identifying unknown clusters and improving confidence in hypothesized cell types.*

2. **Identification of conserved markers for each cluster:** This analysis looks for genes that are differentially expressed/present within each condition first, and then reports those genes that are conserved in the cluster across all conditions. These genes can help to figure out the identity for the cluster. 
	- *Useful with more than one condition to identify cell type markers that are conserved across conditions.*  	
 
3. **Marker identification between specific clusters:** this analysis explores differentially expressed genes between specific clusters. 
	- *Useful for determining differences in gene expression between clusters that appear to be representing the same celltype (i.e with markers that are similar) from the above analyses.*


## Identification of conserved markers in all conditions

Since we have samples representing different conditions in our dataset, our best option is to find conserved markers. Identifying conserved markers allows for identifying only those genes the are significantly differentially expressed relative to the other clusters for all conditions. This function performs differential gene expression testing for a single cluster against all other clusters within each group and then combines the p-values using meta-analysis methods from the MetaDE R package.

Before we start our marker identification we need to explicitly set our default assay. Because of the nature of how `FindConservedMarkers()` works (i.e finding DE within each group and then looking for conservation), we need to use the **original counts and not the integrated data**.

```r
DefaultAssay(combined) <- "RNA"
```

This function is used on multiple samples in lieu of `FindAllMarkers()`. You could run it on all clusters if you wanted to, but it takes a while to run, so we are just going to **run it on the unknown clusters 17 and 20**.

The function `FindConservedMarkers()`, has the following structure:

**`FindConservedMarkers()` syntax:**

```r
FindConservedMarkers(seurat_obj,
                     ident.1 = cluster,
                     grouping.var = "group",
                     only.pos = TRUE)
```

The function **accepts a single cluster at a time**, so if we want to have the function run on all clusters, then we can use the `map` family of functions to iterate across clusters. 

Since these functions will **remove our row names** (gene names), we need to transfer the row names to columns before mapping across clusters. We also need a column specifying **to which cluster the significant genes correspond**.

To do that we will **create our own function** to:

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

[Click here for next lesson]()

***


*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
