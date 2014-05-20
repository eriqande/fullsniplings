

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
#' @examples
#' snp_genos <- get_snp_genos(fs_dev_test_data$plain_snp_data)
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
#' @examples
#' snp_genos <- get_snp_genos(fs_dev_test_data$plain_snp_data)
#' snp_indics <- genos_to_indicators(g = snp_genos$mat)
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
#' @examples 
#' snp_genos <- get_snp_genos(fs_dev_test_data$plain_snp_data)
#' snp_indics <- genos_to_indicators(g = snp_genos$mat)
#' geno_counts <- count_genos(snp_indics)
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
#' @examples 
#' snp_genos <- get_snp_genos(fs_dev_test_data$plain_snp_data)
#' snp_indics <- genos_to_indicators(g = snp_genos$mat)
#' geno_counts <- count_genos(snp_indics)
#' afreqs <- alle_freqs(geno_counts)
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
  if(any(a>1.0000000001 | a<0)) stop("Allele freqs >1 or <0 in a")
  if(any(abs(colSums(a)-1) > 1e-08)) stop("Allele freqs not summing to 1 in a")
  
  ret <- rbind(a[1,] * a[1,], 2 * a[1,] * a[2,], a[2,]*a[2,])
  dimnames(ret) <- list(Genos=1:3, Loci=colnames(a))
  ret
}


#' creates a matrix of genotype probabilities for two unrelated individuals
#'
#' Given a 3 x L array of genotype frequencies, this function returns
#' a 9 x L array of the joint genotype frequencies for two unrelated
#' individuals.  The natural represenation of this would be a 3 x 3 x L
#' array, so the class of the return value is a geno_qty_array
#' @param gf a 3 x L array of genotype frequencies
#' @return a geno_qty array.
#' @export
unrelated_pair_gfreqs <- function(gf) {
  if(any(gf>1.0000000001 | gf<0)) stop("Genotype freqs >1 or <0 in gf")
  if(any(abs(colSums(gf)-1) > 1e-08)) stop("Genotype freqs not summing to 1 in gf")
  if(dim(gf)[1] != 3) stop("gf must be a matrix with three rows")
  ret <- sapply(1:ncol(gf), function(x) outer(gf[,x],gf[,x]))
  dimnames(ret) <- list(UnrelatedPairGenos=paste(0:2, rep(0:2, each=3), sep="-"), Loci=colnames(gf))
  attr(ret, "gqa_dimnames") <- list(Ind1 = 0:2, Ind2 = 0:2, Loci = colnames(gf))
  attr(ret, "gqa_dim") <- c(3, 3, ncol(gf))
  class(ret) <- "geno_qty_array"
  ret
}




#' returns genotype probabilities of a pair of full siblings
#' 
#' Given that the parents are unrelated and have expected genotype
#' frequencies given in the 3 x L array gf, and given genotypin error
#' rates specified for each locus mu (will recycle) this function returns
#' the probablities of the 9 possible genotypes (at all L loci) of
#' a pair of full siblings in the population. This can be done entirely in R
#' but it would be hard to come back to it and understand what it is doing,
#' so I am going to do much of it with a call to an Rcpp function
#' @param gf array of genotype frequencies
#' @param mu genotyping error rates
#' @export
full_sibling_pair_gfreqs <- function(gf, mu=0) {
  pp <- unrelated_pair_gfreqs(gf)  # joint parent pair genotype probs
  tp <- trans_probs()  # transmission probs from parents to a single offpring.  This is 3 x 3 x 3
  ge <- lik_array_from_simple_geno_err(ncol(gf), mu)
  
  # here make a call to the Rcpp function
  ret <- matrix(C_full_sibling_pair_gfreqs(ncol(gf), 3, pp, tp, ge), nrow=9)
  
  dimnames(ret) <- list(FullSibPairGenos=paste(0:2, rep(0:2, each=3), sep="-"), Loci=colnames(gf))
  attr(ret, "gqa_dimnames") <- list(Ind1 = 0:2, Ind2 = 0:2, Loci = colnames(gf))
  attr(ret, "gqa_dim") <- c(3, 3, ncol(gf))
  class(ret) <- "geno_qty_array"
  
  ret
}




#' for each individual, find the other individuals with a sib-pair log likelihood exceeding a certain amount
#' 
#' DEPRECATED.  I keep this around to remind me how slow R is.  It is better
#' to use high_logl_pairs(), implemented in Rcpp, in practice.
#' @param FSP full-sibling pair geno freqs as returned by full_sibling_pair_gfreqs()
#' @param UPF unrelated pair geno freqs as would be returned by unrelated_pair_gfreqs()
#' @param GI  genotypes of all the individuals as a snp_indics array
#' @param loglV  the log-likelihood value that must be exceeded to be included
#' @export
find_high_logl_sib_pairs <- function(FSP, UPF, GI, loglV) {
  # have this while developing
  # FSP <- full_sibling_pair_gfreqs(ret$Gfreqs, mu=0.005)  
  # UPF <- unrelated_pair_gfreqs(gf = ret$Gfreqs)
  # GI <- Vars$snp_indics
  # loglV <- -5.0
  
  L <- dim(GI)[2]
  N <- dim(GI)[3]
  
  # make an array in which each column is the genotype indicators at all loci of one individual'
  # replicated at each locus three times (as if they were individual 1 in FSP)
  v <- matrix(GI, nrow=3)
  g1 <- matrix(rbind(v,v,v), ncol=N)
  
  # now we make an array that is replicated at each element 3 times, as would be individual 2's
  # genotype indicator
  g2 <- matrix(rep(GI, each=3), ncol=N)
  
  # here are the logs of the full sibling probs and the unrelated pair probs
  logFSP <- as.vector(log(FSP))
  logUPF <- as.vector(log(UPF))
  g2 <- as.vector(g2)  # we are only going to use this as a vector, ever.
  
  boing <- lapply(1:N, function(i) {
    gg <- g1[,i] * g2
    pp <- colSums(matrix(gg * logFSP, ncol=N), na.rm=T) -
      colSums(matrix(gg * logUPF, ncol=N), na.rm=T)
    w <- which(pp>loglV)
    w <- w[w != i]  # drop the pair between i and himself!
#    list(idx = w-1, vals = pp[w])  # return the base-0 indices of the fish that are not clearly non-siblings.
    w-1  # I don't need to carry around all the values.
  })
  
  
  
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
      TrueGeno=0:2,
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



#' given some genotype likelihoods and a prior, compute the posterior
#' 
#' @param Lik an array of genotype likelihoods.  These will commonly be "something by N" where
#' N is the number of individuals or marriages.  In that case, the "something" will be all the
#' possible genotype states in such an entity.  (For example 9 x L for a marriage node)
#' @param Pri prior on genotypes.  I will define this vaguely for now.  It must, of course
#' correspond strucurally to Lik.  For example, for marriages, one will expect it to of 
#' length 9 x L.  For individuals likelihoods, 3 x L, an so forth.
#' @param nStates Number of genotypic states at each locus.  So, basically, Pri should be 
#' of length nStates * L, and Lik should be of length nStates * L * N.  I could write 
#' a function that figured out what nStates was from the array dimensions, but I actually
#' think this way is more flexible.
#' @details This could all be done within Rcpp and would probably be faster that way, but I am
#' just prototyping this at the moment.
#' @export
geno_posterior <- function(Lik, Pri, nStates) {
  if(length(Lik) %% length(Pri) != 0) stop("Length of Lik must be a multiple of the length of Pri.")
  
  pp <- as.vector(Pri) * Lik
  normo <- rep(colSums(matrix(pp, nrow=nStates)), each = nStates)
  
  pp/normo
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