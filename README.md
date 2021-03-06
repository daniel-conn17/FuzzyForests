# `fuzzyforest`

`fuzzyforest` is an extension of random forests designed to yield less biased
variable importance rankings when features are correlated with one another.
The algorithm requires that features be partitioned into seperate groups
or modules such that the correlation within groups are large and the 
correlation between groups is small.  `fuzzyforest` allows for easy integration
the package `WGCNA`.

* the latest released version can be downloaded from CRAN: 
  `install.packages("fuzzyforest")`

To enable use of the full functionality of `fuzzyforest` packages `WGCNA`
must be installed.  However, `WGCNA` requires the installation of a few
packages form bioConductor.  To install `WGCNA`, type the following lines
into the console:
```{r}
install.packages("BiocManager")
BiocManager::install("WGCNA")
```

If further issues with the installation of `WGCNA` arise see the `WGCNA`
website: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/


# Acknowledgments

This work is partially supported through NSF grant IIS 1251151.
