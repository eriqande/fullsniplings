
#' S3 class for arrays holding probabilities and likelihoods of different genotypes at different loci
#' 
#' This is the main S3 class for arrays that hold things like the likelihoods of the pair of parental
#' genotypes given all their offspring, or the likelihood of an individual's genotypes given
#' his observed phenotype.  Or, I suppose they could even be indicators of SNP genotypes in an 
#' individual, or normalized probabilities, etc.  Several different objects inherit from
#' this one.  I mostly have this because I wanted to learn more about S3 objects and I wanted
#' to have an easy way of printing these to get information.  
#' 
#' @details Every object that inherits from this class is a numeric (or possibly an integer)
#' matrix.  However, they are typically all objects that would be more naturally represented
#' by higher dimensional arrays, but we don't do that because the operations that we will
#' end up doing on them are typically column and row multiplies, etc, and it is easier in
#' Rcpp to just deal with matrices.  
#' 
#' But, it would be nice to know about the underlying natural higher dimensional array.  These
#' properties are stored in these objects.  The attributes these objects have are:
#' \describe{
#'  \item{dim}{They have a 2-D dim b/c they are matrices}
#'  \item{class}{They may have variable classes, but they all inherit from the \code{geno_qty_array}}
#'  \item{gqa_dim}{The dimension attributes that would be natural to have for the object.}
#'  \item{gqa_dimnames}{The dimnames on those \code{gqa_dim}s.  Note that these are not technically
#'  required, but printing won't work very well without them, and you really should have some
#'  names on those dimensions.}
#'  \item{gqa_description}{A little blurb about what kinds of quantities this object is storing.} 
#' }
#' @name geno_qty_array
#' NULL





#' a method for printing objects inheriting from class geno_qty_array
#' 
#' Note that this method is called just by passing
#' an object of class marriage_geno_lik_array to the print function.
#' @param g An object of class \code{\link{geno_qty_array}}
#' @export
print.geno_qty_array <- function(g) {
  cat("geno_qty_array object of class:  ", class(g), "\n\n")
  cat("Description: ", attr(g, "gqa_description"), "\n\n")
  cat(paste("Matrix representation has dim (", paste(dim(g), collapse=", "), ")", sep=""), "\n")
  cat("Natural dimensional reprentation has dims:\n\n")
  print(sapply(attr(g, "gqa_dimnames"), length) )
  
  cat("\nTo view the object in its natural dimensional representation you can pass\n")
  cat("it to the function gqa_natural() and then subsript that. For example, if g were\n") 
  cat("the name of the object you could do the following to see all of it:\n\n\t")
  cat(paste("gqa_natural(g)[", paste("1:", attr(g, "gqa_dim"), collapse=", ", sep=""), "]", sep=""))
  cat("\n")
}



#' return a geno_qty_array in its natural representation with its dimanames for viewing
#' @param g An object that inherits from the \code{\link{geno_qty_array}} class
#' @export
gqa_natural <- function(g) {
  dim(g) <- attr(g, "gqa_dim")
  dimnames(g) <- attr(g, "gqa_dimnames")
  g
}

