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
```R
# download, build and install the package
devtools::install_github(repo="fullsniplings", username="eriqande")
```

## Running a Simple Data set
Take a look at the first few lines of a SNP data set that has 1157 Chinook typed at 95
SNPs:
```R
head(fs_dev_test_data$chinook_full_sibs_genos)
```
It is simple two column format (each locus gets two alleles) with the rownames being
the IDs of the fish.  Missing data are denote by NA, as is typical in R.

The main function for running the MCMC is `run_mcmc`.  It is not fully documented, yet, but
you can do:
```R
?run_mcmc
```
to learn about it.  The example shows how to do a short run.

If you want to do a longer run (about 5 minutes) you can do
```R
quick_run  <- run_mcmc(fs_dev_test_data$chinook_full_sibs_genos, burn_in = 50, num_sweeps = 200)
```

## Building a vignette that shows how I used fullsniplings for my talk at the Coastwide meeting in June of 2014
To do this, you will need V2 of the `rmarkdown` package.  You can get that  either
by using the latest preview version of Rstudio: http://www.rstudio.com/products/rstudio/download/preview/
or by installing rmarkdown from github like this:
```R
devtools::install_github(repo="rmarkdown", username="studio")
```
Then, once you have done that you should be able to do 
```R
rmarkdown::render(system.file("dev_vignettes/full-sibling-chinook-mcmc.Rmd", package="fullsniplings"), output_format = "html_document", output_file = "~/full-sibling-chinook-mcmc.html")
```
That will put the output in your home directory in the file `full-sibling-chinook-mcmc.html` as long as your system
recognizes the `~`.  If not, just put an absolute path of your choice in for the `output_file` parameter in the 
`render` function.

Note that this doesn't run the MCMC, it just grabs some stored results.  But you can go into the code of 
```R
system.file("dev_vignettes/full-sibling-chinook-mcmc.Rmd", package="fullsniplings")
```
and uncomment the appropriate lines to run it all.
