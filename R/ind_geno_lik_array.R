

#' given genotyping error and observed genotypes, return an array of likelihoods for all the individuals
#' 
#' @param SI  a G x L x N array of snp genotype indicators, of the sort that would come out of 
#'  \code{\link{genos_to_indicators}}.  G is the number of genotypic states, L the number of loci,
#'  and N the number of individuals.
#' @param mu array of per-gene-copy genotyping error rates.  Gets recycled as need be.
#' @export
get_indiv_geno_lik <- function(SI, mu) {
  if(length(dim(SI)) != 3) stop("argument SI must be a three-dimensional array.")
  G <- dim(SI)[1] # number of genotypic states
  L <- dim(SI)[2]
  N <- dim(SI)[3]
  u.mat <- lik_array_from_simple_geno_err(L, mu)
  
  # now we want to pick columns out of u.mat according to the observed genotype
  # of individuals.  We can do this with the indicators in SI, interpreting them as
  # logical vectors. However, every NA will pick out a whole column of u.mat, and we only
  # want to pick out 1 of those, so we turn the top and bottom NA into a 0:
  SI[1,,][is.na(SI[1,,])] <- 0
  SI[3,,][is.na(SI[3,,])] <- 0
  
  # now we replicate u.mat N times and pick out columns of it like so:
  ret <- matrix(rep(u.mat, N), nrow=3)[, as.logical(SI)]
  dim(ret) <- c(G * L, N)  # just make sure this is set
  
  # make it an object of class ind_geno_lik_array which inherits from geno_qty_array
  attr(ret, "gqa_dim") <- dim(SI)
  attr(ret, "gqa_dimnames") <- dimnames(SI)
  class(ret) <- c("ind_geno_lik_array", "geno_qty_array")
  
  ret[is.na(ret)] <- 1 # if no data observed, likelihood of true geno is constant (set to 1) 
  ret
}



