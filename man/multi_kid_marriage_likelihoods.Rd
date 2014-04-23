\name{multi_kid_marriage_likelihoods}
\alias{multi_kid_marriage_likelihoods}
\title{compute the marriage node likelihoods given offspring specified in a full-sibling list}
\usage{
multi_kid_marriage_likelihoods(S, PK)
}
\arguments{
  \item{S}{a list of vectors that give the indices (base 0)
  of the individuals in the full sibling groups.}

  \item{PK}{per-kid marriage likelihoods.  This must be of
  class \code{\link{marriage_geno_lik_array}}, which is
  just a matrix underneath with G x G x L rows and N
  columns.}
}
\description{
compute the marriage node likelihoods given offspring specified in a full-sibling list
}

