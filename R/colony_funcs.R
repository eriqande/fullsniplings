

#### Functions for Reading Colony Output

#' Slurps a colony inferred sibship file into a list 
#' 
#' It also  sorts it by size of sibship and 
#' probability.
#' @param file the name of the colony output file.  It should be one of the ones with 
#' the extension \code{.BestFSFamily}. 
#' @param skip a vector of lines in the file that
#' should be skipped or ignored. (like maybe the first line...)
#' @return A list of two named componenets:
#' \describe{
#'  \item{probs}{Colony sibship approximated posterior probababilities}
#'  \item{sibs}{Vectors of the individuals in each sibship}
#' }
#' @export
SlurpSibships <- function(file, skip=integer(0)) {
  csibs <- readLines(con=file)
  if(length(skip)>0) {
    csibs <- strsplit(csibs[-skip], split=" +")
  } else {
    csibs <- strsplit(csibs, split=" +")
  }
  cs.probs <- sapply(csibs, function(x) x[3])
  cs.sibs <- lapply(csibs, function(x) x[-(1:3)])
  ord <- order( sapply(cs.sibs, length), cs.probs, decreasing=c(T,F))
  list(probs=cs.probs[ord], sibs=cs.sibs[ord])
}


#' adds the names of the members of the sibset if available 
#' 
#' Also adds the indices in the data set of those members, if they are 
#' available.  
#' @param S The output of SlurpSibships.
#' @param N  a vector whose names attribute is the names in S
#' and whose values are whatever names you want for the sibs.
#' @param D  a vector of indices with names being those in S
#' @export
NamizeSibs <- function(S, N=NULL, D=NULL) {
  if(is.null(N[1])) {
    nn<-NULL
  }
  else {
    nn<-lapply(S$sibs, function(x) {r<-N[x]; names(r)<-NULL; r})
    names(nn)<-NULL
  }
  if(is.null(D[1])) {
    dd<-NULL
  }
  else {
    dd<-lapply(S$sibs, function(x) {r<-D[x]; names(r)<-NULL; r})
    names(dd)<-NULL
  }
  c(S, list(sibs.names=nn, sibs.idx=dd))
}



