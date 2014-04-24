\name{full_sibling_pair_greqs}
\alias{full_sibling_pair_greqs}
\title{returns genotype probabilities of a pair of full siblings}
\usage{
full_sibling_pair_greqs(gf, mu = NULL)
}
\description{
Given that the parents are unrelated and have expected genotype
frequencies given in the 3 x L array gf, and given genotypin error
rates specified for each locus mu (will recycle) this function returns
the probablities of the 9 possible genotypes (at all L loci) of
a pair of full siblings in the population
}

