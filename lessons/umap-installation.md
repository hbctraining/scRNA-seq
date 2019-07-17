### Installing and using UMAP


To visualize the cell clusters, there are a few different dimensionality reduction techniques that can be helpful. The most popular methods include t-distributed stochastic neighbor embedding (t-SNE) and **Uniform Manifold Approximation and Projection (UMAP)** techniques.

In the Seurat package there is a function to use the UMAP visualization (`RunUMAP()`), however it does require the user to  **first install the `umap-learn` python package.**  There are a few ways to go about this, and we have them listed below in the order in which you should try them.

1. Since this is a Python package, you will first need Python installed on your laptop (and the correct version). **Download and install Miniconda** using the [link provided](https://docs.conda.io/en/latest/miniconda.html). Be sure to choose the installer for Python 3.7 and the appropriate OS you have on your laptop.


2. **In the R console**, type in the following:

```r
system("which python")
```

You should see a path output that identifies where the executable for Python is stored, for example:

```r
/Users/hbctraining/miniconda3/bin/python

```

**If you see this please continue to #4.


3. If you **don't see miniconda in that path**, you will want to set the path to use your recent installation of Python.

Miniconda will have been installed in your home directory. To find the full path to your home directory, type out:

```r
system("echo $HOME")

> /Users/hbctraining/

```

Then concatenate your home path with `/miniconda3/bin/python` and set this in your environment using the `use_python()` function from the `reticulate` library:

```r
library(reticulate)
use_python(python = "/Users/hbctraining/miniconda3/bin/python", required = T)

```

4. Now you can install the `umap-learn` package using teh `py_install()` function, also from your `reticulate` library:

```r
library(reticulate) # if not already loaded
py_install("umap-learn")

```

If you restart your R session this will now take effect and the `RunUMAP()` *should* work without errors.

*Keep an eye out for any error messages!!!* 

