

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


#' creates a vector of genotypes likelihoods from the allele freqs
#' 
#' Given a 2 x L array of allele frequencies, this function returns
#' a 3 x L array of genotype frequencies expected under Hardy-Weinberg
#' proportions.
#' @param a the 2 x L array of allele frequencies.  
#' @return a 3 x L array of genotypes frequencies
#' @export
gfreqs_from_afreqs <- function(a) {
  if(any(a>1 | a<0)) stop("Allele freqs >1 or <0 in a")
  if(any(colSums(a)>1)) stop("Allele freqs summing to more than 1 in a")
  
  ret <- rbind(a[1,] * a[1,], 2 * a[1,] * a[2,], a[2,]*a[2,])
  dimnames(ret) <- list(Genos=1:3, Loci=colnames(a))
  ret
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


#' a function to prep all the variables up through the marriage likelihoods
#' 
#' @param DF the data frame of SNP data
#' @param geno_error_rates  vector of genotyping error rates for each locus (recycled if necessary)
#' @export
prep_all_variables <- function(genos, geno_error_rates) {
  snp_genos <- get_snp_genos(genos)
  snp_indics <- genos_to_indicators(g = snp_genos$mat)
  geno_counts <- count_genos(snp_indics)
  afreqs <- alle_freqs(geno_counts)
  geno_liks <- get_indiv_geno_lik(SI = snp_indics, mu = geno_error_rates)
  pk_marriage_liks <- per_kid_marriage_likelihoods(geno_liks, trans_probs())
  
  list(
    snp_genos = snp_genos,
    snp_indics = snp_indics,
    geno_counts = geno_counts,
    afreqs = afreqs,
    geno_liks = geno_liks,
    pk_marriage_liks = pk_marriage_liks
    )
}