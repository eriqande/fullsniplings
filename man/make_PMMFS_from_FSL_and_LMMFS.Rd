\name{make_PMMFS_from_FSL_and_LMMFS}
\alias{make_PMMFS_from_FSL_and_LMMFS}
\title{create a "Posterior Matrix of Marriages given Full Siblings" from a full sibling
list and LMMFS}
\usage{
make_PMMFS_from_FSL_and_LMMFS(FSL, LMMFS, af)
}
\arguments{
  \item{FSL}{a full sibling list}

  \item{LMMFS}{a Likelihood Matrix of Marriages given Full
  Siblings}

  \item{af}{the vector of allele frequencies}
}
\description{
This is mostly a wrapper for update_marriage_posteriors_in_place.  It assumes that the
genotypes of the parents of each marriage are just drawn from the allele frequencies
under Hardy-Weinberg equilibrium
}

