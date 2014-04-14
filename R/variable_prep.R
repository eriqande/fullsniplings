

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