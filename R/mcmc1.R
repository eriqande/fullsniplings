

#' do a whole MCMC analysis starting from a data set and genotyping error rates
#' 
#' @param genos  The data frame of SNP genotypes
#' @param mu  Genotyping error rates per locus (recycles as necessary)
#' @export
full_sib_mcmc <- function(genos, mu) {
  
  # first get all the data. Store it in a list that we will add other Variables to, as well
  Vars <- prep_all_variables(genos = genos, geno_error_rates = mu)
  
  # get the full sibling list
  FSL <- initialize_sib_list_to_singletons( ncol(Vars$geno_liks) )
  
  # from that list, make the "Individuals' full sibling groups" vector
  IFS <- make_IFS_from_FSL(FSL, ncol(Vars$geno_liks))
  
  # set the LMMI to be Vars$pk_marriage_liks (just to keep notation consistent with my notebook notes)
  LMMI <- Vars$pk_marriage_liks
  
  # get the initial condition for the LMMFS
  LMMFS <- make_LMMFS_from_FSL_and_LMMI(FSL, LMMI)
  
  
}



#' initialize the list of full siblings to a bunch of singletons
#' 
#' @param N number of individuals
#' @export
initialize_sib_list_to_singletons <- function(N) {
  lapply(1:N, function(x) list(LMMI_Idx = x-1, Indivs = x-1))
}


#' create the "individuals full sibships" vector from the current state of the full sibling list
#' 
#' @param FSL the full sibling list
#' @param N the number of individuals that you should have entries for
#' @export
make_IFS_from_FSL <- function(FSL, N) {
  tmp <- do.call(rbind, lapply(FSL, function(x) cbind(x$Indivs, x$LMMI_Idx)))
  rownames(tmp) <- as.character(tmp[,1])
  tmp <- unname(tmp[as.character(0:(N-1)),2])
  if(any(is.na(tmp))) stop("Missing some indices.  Aborting.")
  tmp
}


#' create the "Likelihood Matrix of Marriages given Full Siblings" from a full sibling list
#' and the Likelihood Matrix of Marriages given Individuals
#' 
#' This is mostly a wrapper for update_marriage_likelihoods_in_place
#' @param FSL a full sibling list
#' @param LMMI a Likelihood Matrix of Marriages given Individuals
#' @export
make_LMMFS_from_FSL_and_LMMI <- function(FSL, LMMI) {
  ret <- LMMI  # make it big
  ret[] <- 0.0  # set all values to 0.0, initially. Note that this also makes ret truly separate from LMMI
  
  update_marriage_likelihoods_in_place(FSL, LMMI, ret, 0:(length(FSL)-1))
  ret
}