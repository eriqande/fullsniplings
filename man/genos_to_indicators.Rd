\name{genos_to_indicators}
\alias{genos_to_indicators}
\title{Take an L x N matrix of SNP genos and return a 3 x L x N array of indicators}
\usage{
genos_to_indicators(g)
}
\arguments{
  \item{g}{an L x N matrix of SNP genotypes that are 0, 1,
  2, or NA.}
}
\description{
You actually pass this a snp_genos object and it deals with putting names and
stuff on there.
}

