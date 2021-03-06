% Generated by roxygen2 (4.0.1): do not edit by hand
\name{make_KidProngs_from_FSL_and_PMMFS}
\alias{make_KidProngs_from_FSL_and_PMMFS}
\title{create a KidProngs matrix (matrix of posterior predictives for the next fullsibling in a full sibship)
from a PMMFS and a FSL.}
\usage{
make_KidProngs_from_FSL_and_PMMFS(FSL, PMMFS)
}
\arguments{
\item{FSL}{a full sibling list}

\item{PMMFS}{a Posterior Matrix of Marriages given Full Siblings}
}
\description{
This is mostly a wrapper for update_marriage_node_kid_prongs_in_place.
}

