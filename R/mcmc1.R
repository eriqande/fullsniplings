

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
  
  # get the initial condition for the posteriors on the parent genotypes from each full sibship in LMMFS
  PMMFS <- make_PMMFS_from_FSL_and_LMMFS(FSL, LMMFS, Vars$afreqs)
  
  # now we set the initial conditions on the Kid-Prongs!
  KidProngs <- make_KidProngs_from_FSL_and_PMMFS(FSL, PMMFS)
  
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
#' @param allocateFromLMMI  should the dimension of LMMI be used to allocate memory 
#' for LMMFS.  If FALSE then it won't overallocate.
#' @export
make_LMMFS_from_FSL_and_LMMI <- function(FSL, LMMI, allocateFromLMMI = TRUE) {
  if(allocateFromLMMI == TRUE) {
    ret <- LMMI  # make it big
    ret[] <- 0.0  # set all values to 0.0, initially. Note that this also makes ret truly separate from LMMI
  } else {
    ret <- matrix(0, nrow=nrow(LMMI), ncol=length(FSL))
  }
  
  update_marriage_likelihoods_in_place(FSL, LMMI, ret, 0:(length(FSL)-1))
  
  ret
}



#' create a "Posterior Matrix of Marriages given Full Siblings" from a full sibling
#' list and LMMFS
#' 
#' This is mostly a wrapper for update_marriage_posteriors_in_place.  It assumes that the 
#' genotypes of the parents of each marriage are just drawn from the allele frequencies
#' under Hardy-Weinberg equilibrium
#' @param FSL a full sibling list
#' @param LMMFS a Likelihood Matrix of Marriages given Full Siblings
#' @param af the vector of allele frequencies
#' @export
make_PMMFS_from_FSL_and_LMMFS <- function(FSL, LMMFS, af) {
  ret <- LMMFS
  ret[] <- 0.0
  
  UPG <- unrelated_pair_gfreqs(gfreqs_from_afreqs(af))  # this is our prior on parent genotypes in the simplest case

  update_marriage_posteriors_in_place(
        FSL, 
        LMMFS, 
        ret, 
        UPG, 
        9, 
        0:(length(FSL)-1)
  )
  ret
}


#' create a KidProngs matrix (matrix of posterior predictives for the next fullsibling in a full sibship) 
#' from a PMMFS and a FSL.  
#' 
#' This is mostly a wrapper for update_marriage_node_kid_prongs_in_place.  
#' @param FSL a full sibling list
#' @param PMMFS a Posterior Matrix of Marriages given Full Siblings
#' @export
make_KidProngs_from_FSL_and_PMMFS <- function(FSL, PMMFS) {
  nL <- nrow(PMMFS) / 9  # number of loci
  nN <- ncol(PMMFS)      # number of individuals
  
  ret <- matrix(0.0, nrow = 3 * nL, ncol = nN)  # make the matrix to return
  
  # then modify ret in place
  update_marriage_node_kid_prongs_in_place(
    FSL, PMMFS, ret, 
    9, 3, trans_probs(), 
    0:(length(FSL)-1))
  
  # and return it.  I really should make it a geno_qty_array, but I am tired of that at the moment
  ret
}