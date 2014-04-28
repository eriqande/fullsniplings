\name{update_marriage_likelihoods_in_place}
\alias{update_marriage_likelihoods_in_place}
\title{Update one row of a marriage likelihood matrix IN PLACE!}
\usage{
update_marriage_likelihoods_in_place(S, PK, ML, bz_idx)
}
\arguments{
  \item{S}{a list of vectors that give the indices (base 0)
  of the individuals in the full sibling groups.}

  \item{PK}{per-kid marriage likelihoods.  This must be of
  class \code{\link{marriage_geno_lik_array}}, which is
  just a matrix underneath with G x G x L rows and N
  columns.}

  \item{ML}{the marriage likelihoods matrix to be modified}

  \item{bz_idx}{An integer vector holding the BASE-0
  indices of the rows of ML to be updated. They will be
  updated according to the current contents of S. //'}
}
\value{
This doesn't return anything.  It modifies ML in place via call be reference.  Our
goal here is to make updates without copying a lot of memory.
}
\description{
The intended use of this is to update the marriage likelihoods after an individual has
been moved from one sibship to another.  So, for example, if you moved an individual
from sibship 40 to sibship 291 (as subscripted in R) then you would want to run this
with \code{bz_idx} equal to \code{c(39, 290)}.
}

