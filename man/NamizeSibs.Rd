% Generated by roxygen2 (4.0.1): do not edit by hand
\name{NamizeSibs}
\alias{NamizeSibs}
\title{adds the names of the members of the sibset if available}
\usage{
NamizeSibs(S, N = NULL, D = NULL)
}
\arguments{
\item{S}{The output of SlurpSibships.}

\item{N}{a vector whose names attribute is the names in S
and whose values are whatever names you want for the sibs.}

\item{D}{a vector of indices with names being those in S}
}
\description{
Also adds the indices in the data set of those members, if they are
available.
}

