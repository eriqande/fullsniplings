# functions for assessing the accuracy of inferred sibships
# but nothing in here yet for computing the partition distance


#' computes the average number of non-sibs in inferred sibships
#' 
#' And also computes the average number of sibs missing.  
#' For each individual i in an inferred sibship this function counts the
#' average number of "errors of commission" (number of individuals 
#' incorrectly declared to be i's sibling.  And it averages over i.
#' It also computes for each individual i the number of his true sibs
#' missing from the sibroup he is inferred to be in.
#' This function operates on a list of sibships
#' @param S A list of inferred sibships where individuals are named by their base-0 indices.
#' @param sg.list A list holding the true siblings of each individual.  If 
#' if \code{Ind} is the base-0 index of an individual then 
#' \code{sg.list[[Ind+1]]} is a vector of the base-0 indices of all of its full siblings
#' @export
ave_err_commission_and_ommission <- function(S, sg.list) {
  
  # define a quick function that computes the mean on a single inferred sibgroup with 
  # elements x
  om_mean <- function(x) {mean(sapply(x, function(z) length(setdiff(c(z, sg.list[[z+1]]), x))))}
  com_mean <- function(x) {mean(sapply(x, function(z) length(setdiff(x, c(z, sg.list[[z+1]])))))}
  
  # then lapply that function over S
  list(AveNumSibsMissing = sapply(S, om_mean), AveNumNonSibs = sapply(S, com_mean))
  
}