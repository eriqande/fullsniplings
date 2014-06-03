# fullsniplings

An R package to infer full-sibling relationships from genetic data.
Geared toward datasets with many individuals and about 100 SNPs which are
treated as if they are unlinked. (While that is clearly false, it might be
reasonable).  The method is based on computing the posterior predictive
distribution of genotypes conditional on currently-inferred sibling groups. 
I will eventually generalize
that to more general pedigree inference, but for now just want to implement
the easy full-sibling case as a prototypical exercise.

## Installing / Building
On a Mac you need to have the developer command line tools, and on a PC you will 
want the RTools to build this package.  You need to have the `Rcpp` package as
well, and also the `devtools` package.  You can get both of those from CRAN.

Here is what you would do:
```{r}
library(devtools)  # load this package for the install_github() function

```