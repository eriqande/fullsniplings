\name{make_LMMFS_from_FSL_and_LMMI}
\alias{make_LMMFS_from_FSL_and_LMMI}
\title{create the "Likelihood Matrix of Marriages given Full Siblings" from a full sibling list
and the Likelihood Matrix of Marriages given Individuals}
\usage{
make_LMMFS_from_FSL_and_LMMI(FSL, LMMI)
}
\arguments{
  \item{FSL}{a full sibling list}

  \item{LMMI}{a Likelihood Matrix of Marriages given
  Individuals}
}
\description{
This is mostly a wrapper for update_marriage_likelihoods_in_place
}

