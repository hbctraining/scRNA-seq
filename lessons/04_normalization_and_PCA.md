---
title: "Single-cell RNA-seq: Theory of normalization and PCA"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Monday, July 15th, 2019
---

Approximate time: 90 minutes

## Learning Objectives 

* Normalizing counts for accurate comparison between cells
* Evaluating similarity in cellular gene expression using Principal Components Analysis

After attaining our high quality single cells, the next step in the single-cell RNA-seq (scRNA-seq) analysis workflow is to perform clustering. The goal of clustering is to separate different cell types into unique clusters of cells. To perform clustering, we determine the genes that are most different in their expression between cells. Then, we use these genes to determine which correlated genes sets are responsible for the largest differences in expression between cells.


<p align="center">
<img src="../img/sc_workflow.png" width="800">
</p>


### Count normalization

Count normalization is necessary to make accurate comparisons of gene expression between cells. The counts of mapped reads for each gene is proportional to the expression of RNA ("interesting") in addition to many other factors ("uninteresting"). Normalization is the process of scaling raw count values to account for the "uninteresting" factors. In this way the expression levels are more comparable between and/or within cells.

The main factors often considered during normalization are:
 
 - **Sequencing depth:** Accounting for sequencing depth is necessary for comparison of gene expression between cells. In the example below, each gene appears to have doubled in expression in cell 2, however this is a consequence of cell 2 having double the sequencing depth.
 
    <img src="../img/sequencing_depth.png" width="400">
    
    Each cell in scRNA-seq will have a differing number of reads associated with it. Some cells may have double the reads of other cells, so to accurately compare expression between cells, it is necessary to normalize for sequencing depth.
 
 - **Gene length:** Accounting for gene length is necessary for comparing expression between different genes within the same cell. The number of reads mapped to a longer gene can appear to have equal count/expression as a shorter gene that is more highly expressed. 
 
    <img src="../img/length_of_gene.png" width="400">
 
	In scRNA-seq analysis, we will be comparing the expression of different genes within the cells to cluster the cells. If using a 3' or 5' droplet-based method, the length of the gene will not affect the analysis because only the 5' or 3' end of the transcript is sequenced. However, if using full-length sequencing, the transcript length should be accounted for.

### Principal Component Analysis (PCA)

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). Details regarding PCA are given below (based on [materials from StatQuest](https://www.youtube.com/watch?v=_UVHneBUBW0), and if you would like a more thorough description, we encourage you to explore [StatQuest's video](https://www.youtube.com/watch?v=_UVHneBUBW0). 

If you had two samples and wanted to plot the counts of one sample versus another, you could plot the counts of one sample on the x-axis and the other sample on the y-axis as shown below:

<img src="../img/PCA_2sample_genes.png" width="600">

You could draw a line through the data in the direction representing the most variation, which is on the diagonal in this example. The maximum variation in the data is between the two endpoints of this line.  

We also see the genes vary somewhat above and below the line. We could draw another line through the data representing the second most amount of variation in the data. 


<img src="../img/PCA_2sample_variation1.png" width="600">

The genes near the ends of the line, which would include those genes with the highest variation between samples (high expression in one sample and low expression in the other), have the greatest influence on the direction of the line. 

<img src="../img/PCA_2sample_variation2.png" width="600">

For example, a small change in the value of *Gene C* would greatly change the direction of the line, whereas a small change in *Gene A* or *Gene D* would have little affect.

<img src="../img/PCA_2sample_variation3.png" width="900">

We could just rotate the entire plot and view the lines representing the variation as left-to-right and up-and-down. We see most of the variation in the data is left-to-right; this is and the second most variation in the data is up-and-down. These axes that represent the variation are "Principal Components", with PC1 representing the most variation in the data and PC2 representing the second most variation in the data. 

If we had three samples, then we would have an extra direction in which we could have variation. Therefore, if we have *N* samples we would have *N*-directions of variation or principal components.

<img src="../img/PCA_2sample_rotate.png" width="300">

We could give quantitative scores to genes based on how much they influence PC1 and PC2. Genes with little influence would get scores near zero, while genes with more influence would receive larger scores. Genes on opposite ends of the lines have a large influence, so they would receive large scores, but with opposite signs.

<img src="../img/PCA_2sample_influence.png" width="600">

To generate a score per sample, we combine the read counts for all genes. To calculate the scores, we do the following:
	
	Sample1 PC1 score = (read count * influence) + ... for all genes
	
Using the counts in the table for each gene (assuming we had only 4 genes total) we could calculate PC1 and PC2 values for each sample as follows:

	Sample1 PC1 score = (4 * -2) + (1 * -10) + (8 * 8) + (5 * 1) = 51
	Sample1 PC2 score = (4 * 0.5) + (1 * 1) + (8 * -5) + (5 * 6) = -7
	
	Sample2 PC1 score = (5 * -2) + (4 * -10) + (8 * 8) + (7 * 1) = 21
	Sample2 PC2 score = (5 * 0.5) + (4 * 1) + (8 * -5) + (7 * 6) = 8.5
	
The scores would then be plotted to examine whether the samples exhibit similar variation across all genes:

<img src="../img/PCA_samples.png" width="600">

Since genes with the greatest variation between samples will have the greatest influence on the principal components, we hope our condition of interest explains this variation (e.g. high counts in one condition and low counts in the other). With PC1 representing the most variation in the data and PC2 representing the second most variation in the data, we can visualize how similar the variation of genes is between samples.

The same analysis could be performed with 2 cells, or in our case, thousands of cells.
