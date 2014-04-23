\name{get_snp_genos}
\alias{get_snp_genos}
\title{takes a data frame of SNP alleles and converts them to a matrix of 0, 1, 2, and NA}
\usage{
get_snp_genos(GD)
}
\arguments{
  \item{snp_frame}{The data frame of SNPs}
}
\value{
This returns a list with the following components:
\describe{
   \item{mat}{An L x N matrix of 0, 1, 2, or NA's}
   \item{alle_names}{The names of the alleles.  The first of these is what the "0" allele
         is while the second is the name of the "1" allele.}
 }
}
\description{
Let L be the number of loci and N the number of alleles.  The data frame should
have 2L columns and N rows.  The rownames should be the individual names.  The
locus names should be the odd numbered column names.  The even numbered column names
are not used.
}
\details{
A SNP allele can be denoted by any string or number or NA.  Each SNP must have
no more than two alleles.
}

