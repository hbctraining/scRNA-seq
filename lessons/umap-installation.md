### Installing and using UMAP


To visualize the cell clusters, there are a few different dimensionality reduction techniques that can be helpful. The most popular methods include t-distributed stochastic neighbor embedding (t-SNE) and **Uniform Manifold Approximation and Projection (UMAP)** techniques.

In the Seurat package there is a function to use the UMAP visualization (`RunUMAP()`), however it does require the user to  **first install the `umap-learn` python package.** 

We **first** recommend trying this installation with in R:

```r
library(reticulate)
py_install("umap-learn")

```

If you restart your R session this will now take effect and the `RunUMAP()` will work without errors.

*Keep an eye out for any error and warning messages.* If it did not work for you, you will need to try installing using a Python package management system (i.e one fo the methods below):

**Check if you have Python on your computer**. In the terminal (*not inside the R console*) type out:

```bash
which python
```

If this returns to you a path, you will know that you have it installed. Next **check what version**:

```bash

python --version

```

In order to install `umap-learn` you should be using Python version 3.0 or higher. If this is not the case please let us know.

If you have the correct version of Python installed you can install `umap-learn` using Python PIP. Note, if you have Python version 3.4 or later, PIP is included by default. If not, you will need to download and install PIP [using this link](https://pypi.org/project/pip/).

Once you have PIP installed you can run the command below in your terminal (*not in R*). 

```bash

pip install umap-learn
```

Once it completes, you can go back to RStudio and restart your R session for it to take effect.


> **NOTE:** Version of Python and PIP will matter. If you have Python 3 and you find it is still not working, it might be that your PIP needs an upgrade. If so, try the commands below:
>
> **Mac Users**
> `pip install -U pip`
>
> **Windows Users**
> `python -m pip install -U pip`
> 
