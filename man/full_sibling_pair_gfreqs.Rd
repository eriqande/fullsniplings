\name{full_sibling_pair_gfreqs}
\alias{full_sibling_pair_gfreqs}
\title{returns genotype probabilities of a pair of full siblings}
\usage{
full_sibling_pair_gfreqs(gf, mu = 0)
}
\arguments{
  \item{gf}{array of genotype frequencies}

  \item{mu}{genotyping error rates}
}
\description{
Given that the parents are unrelated and have expected genotype
frequencies given in the 3 x L array gf, and given genotypin error
rates specified for each locus mu (will recycle) this function returns
the probablities of the 9 possible genotypes (at all L loci) of
a pair of full siblings in the population. This can be done entirely in R
but it would be hard to come back to it and understand what it is doing,
so I am going to do much of it with a call to an Rcpp function
}

