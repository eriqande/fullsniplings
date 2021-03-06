% Generated by roxygen2 (4.0.1): do not edit by hand
\name{geno_qty_array}
\alias{geno_qty_array}
\title{S3 class for arrays holding probabilities and likelihoods of different genotypes at different loci}
\description{
S3 class for arrays holding probabilities and likelihoods of different genotypes at different loci
}
\details{
This is the main S3 class for arrays that hold things like the likelihoods of the pair of parental
genotypes given all their offspring, or the likelihood of an individual's genotypes given
his observed phenotype.  Or, I suppose they could even be indicators of SNP genotypes in an
individual, or normalized probabilities, etc.  Several different objects inherit from
this one.  I mostly have this because I wanted to learn more about S3 objects and I wanted
to have an easy way of printing these to get information.

Every object that inherits from this class is a numeric (or possibly an integer)
matrix.  However, they are typically all objects that would be more naturally represented
by higher dimensional arrays, but we don't do that because the operations that we will
end up doing on them are typically column and row multiplies, etc, and it is easier in
Rcpp to just deal with matrices.

But, it would be nice to know about the underlying natural higher dimensional array.  These
properties are stored in these objects.  The attributes these objects have are:
\describe{
 \item{dim}{They have a 2-D dim b/c they are matrices}
 \item{class}{They may have variable classes, but they all inherit from the \code{geno_qty_array}}
 \item{gqa_dim}{The dimension attributes that would be natural to have for the object.}
 \item{gqa_dimnames}{The dimnames on those \code{gqa_dim}s.  Note that these are not technically
 required, but printing won't work very well without them, and you really should have some
 names on those dimensions.}
 \item{gqa_description}{A little blurb about what kinds of quantities this object is storing.}
}
}

