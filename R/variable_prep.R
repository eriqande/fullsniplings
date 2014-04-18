

#' takes a data frame of SNP alleles and converts them to a matrix of 0, 1, 2, and NA
#' 
#' Let L be the number of loci and N the number of alleles.  The data frame should 
#' have 2L columns and N rows.  The rownames should be the individual names.  The
#' locus names should be the odd numbered column names.  The even numbered column names
#' are not used.
#' 
#' A SNP allele can be denoted by any string or number or NA.  Each SNP must have
#' no more than two alleles.
#' @param snp_frame  The data frame of SNPs
#' @return This returns a list with the following components:
#' \describe{
#'    \item{mat}{An L x N matrix of 0, 1, 2, or NA's}
#'    \item{alle_names}{The names of the alleles.  The first of these is what the "0" allele
#'          is while the second is the name of the "1" allele.}
#'  }
#' @export
get_snp_genos <- function(GD) {
  if(ncol(GD) %% 2 != 0) stop("GD must have an even number of columns")
  L <- ncol(GD) / 2
  
  get_alles <- function(x) {
    levs <- levels(as.factor(c(GD[,2*x-1], GD[,2*x])))
    a <- factor(GD[,2*x-1], levels=levs)
    b <- factor(GD[,2*x], levels=levs)
    if(length(levels(a)) > 2) stop(paste("Too many alleles at locus", x))
    list(alles = levels(a), vals = as.integer(a)-1 + as.integer(b)-1)
  }
  
  tmp <- lapply(1:L, get_alles)
  mat <- do.call(rbind, lapply(tmp, function(x) x$vals))
  alles <- lapply(tmp, function(x) x$alles)
  names(alles) <- colnames(GD)[c(T,F)]
  rownames(mat) <- colnames(GD)[c(T,F)]
  colnames(mat) <- rownames(GD)
  
  list(mat = mat, alle_names = alles)
}




#' Take an L x N matrix of SNP genos and return a 3 x L x N array of indicators
#' 
#' You actually pass this a snp_genos object and it deals with putting names and 
#' stuff on there.
#' 
#' @param g  an L x N matrix of SNP genotypes that are 0, 1, 2, or NA.
#' @export
#' 
genos_to_indicators <- function(g) {
  
  x <- rbind( c(1,0,0),
              c(0,1,0),
              c(0,0,1))
  
  y <- x[, g+1]
  dim(y) <- c(3, nrow(g), ncol(g))
  dimnames(y) <- list(Genos = 0:2, 
                      Loci = rownames(g),
                      Indivs = colnames(g))
  y
}



#' counts up genotypes from a snp_indicators array
#' 
#' @param x a snp indicators array (3 x L x N)
#' @export
count_genos <- function(S) {
  d <- dim(S)
  x <- S
  dim(x) <- c(d[1] * d[2], d[3])
  y <- matrix(rowSums(x, na.rm=T), nrow=3)
  dimnames(y) <- dimnames(S)[1:2]
  y
}


#' counts alleles from a geno_counts matrix
#' 
#' @param x the geno counts matrix
#' @param proportion  if false it returns the counts. if true the proportions (summing to one)
#' @param smidge small amount to add to count of each allele (could be useful if there are zero counts on some alleles)
#' @return A matrix with 2 rows (corresponding to alleles 0 and 1) and L columns.
#' @export
alle_freqs <- function(x, proportion=T, smidge=0.5) {
  y <- rbind(colSums(x * c(2, 1, 0) , na.rm=T), colSums(x * c(0, 1, 2), na.rm=T) )
  dimnames(y)[[1]] <- 0:1
  names(dimnames(y)) <- c("Alleles", "Loci")
  if(proportion==TRUE) {
    return((y+smidge) / rbind(colSums(y+smidge), colSums(y+smidge)))
  }
  y
}



#' makes L 3x3 tables of genotype likelihoods
#' 
#' @param L the desired number of loci
#' @param mu vector of per-gene-copy genotyping error rates, one per locus.  Will recycle (or be truncated) to be of length L.
#' @export
lik_array_from_simple_geno_err <- function(L, mu) {
  mu <- rep(mu, length.out=L)  # recycle as need be
  y <- matrix(sapply(mu, function(u) {
    c(
      (1-u)^2,  u*(1-u),  u^2,
      2*u*(1-u), u^2 + (1-u)^2, 2*u*(1-u),
      u^2,  u*(1-u),  (1-u)^2
      )
  }), nrow=3)
  
  dimnames(y) <- list(
      TrueGeno=1:3,
      ObsGeno=paste("mu_", rep(1:L, each=3), ".", 0:2, sep="")
    )
  y
}


#' given genotyping error and observed genotypes, return an array of likelihoods for all the individuals
#' 
#' @param SI  a 3 x L x N array of snp genotype indicators, of the sort that would come out of 
#'  \code{\link{genos_to_indicators}}.
#' @param mu array of per-gene-copy genotyping error rates.  Gets recycled as need be.
#' @export
get_indiv_geno_lik <- function(SI, mu) {
  if(length(dim(SI)) != 3) stop("argument SI must be a three-dimensional array.")
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
  dim(ret) <- dim(SI)
  dimnames(ret) <- dimnames(SI)
  ret[is.na(ret)] <- 1 # if no data observed, likelihood of true geno is constant (set to 1) 
  ret
}



#' returns a matrix of the 27 transmission probs between a pair of parents and an offspring
#' 
#' There are 9 possible genotype configurations among the parents and 3 in the kid, this
#' function just enumerates those and puts the probabilities into a 3 x 3 x 3 matrix and
#' returns it. It initializes all values to 0, then only modifies the non-zero ones.
#' @export
trans_probs <- function() {
  ret <- rep(0.0,27)
  dim(ret) <- c(3,3,3)
  dimnames(ret) <- list(ParentOne = paste(0:2), ParentTwo = paste(0:2), Kid = paste(0:2))
  
  ret["0", "0", "0"] = 1.0
  ret["0", "1", "0"] = 0.5
  ret["0", "1", "1"] = 0.5
  ret["0", "2", "1"] = 1.0
  ret["1", "0", "0"] = 0.5
  ret["1", "0", "1"] = 0.5
  ret["1", "1", "0"] = 0.25
  ret["1", "1", "1"] = 0.5
  ret["1", "1", "2"] = 0.25
  ret["1", "2", "1"] = 0.5
  ret["1", "2", "2"] = 0.5
  ret["2", "0", "1"] = 1.0
  ret["2", "1", "1"] = 0.5
  ret["2", "1", "2"] = 0.5
  ret["2", "2", "2"] = 1.0
  
  ret
}