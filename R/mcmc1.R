

#' initilize all the variables to do a full sibling analysis given genos and mu
#' 
#' This returns everything in a named list suitable for do.calling 
#' the MCMC sweep function.
#' @param genos  The data frame of SNP genotypes
#' @param mu  Genotyping error rates per locus (recycles as necessary)
#' @return This returns a list with components that are documented as the parameters that 
#' get passed to \link{\code{gibbs_update_one_indiv_in_place}}
#' @export
full_sib_mcmc_initialize <- function(genos, mu) {
  
  message("Starting mcmc variable initialization at ", date())
  
  # first get all the data. Store it in a list that we will add other Variables to, as well
  Vars <- prep_all_variables(genos = genos, geno_error_rates = mu)
  
  ret <- list()
  
  # get the full sibling list
  ret$FSL <- initialize_sib_list_to_singletons( ncol(Vars$geno_liks) )
  
  # from that list, make the "Individuals' full sibling groups" vector
  ret$IFS <- as.integer(make_IFS_from_FSL(ret$FSL, ncol(Vars$geno_liks)))
  
  # set the LMMI to be Vars$pk_marriage_liks (just to keep notation consistent with my notebook notes)
  ret$LMMI <- Vars$pk_marriage_liks
  
  # get the initial condition for the LMMFS
  ret$LMMFS <- make_LMMFS_from_FSL_and_LMMI(ret$FSL, ret$LMMI)
  
  # get the initial condition for the posteriors on the parent genotypes from each full sibship in LMMFS
  ret$PMMFS <- make_PMMFS_from_FSL_and_LMMFS(ret$FSL, ret$LMMFS, Vars$afreqs)
  
  # now we set the initial conditions on the Kid-Prongs!
  ret$KidProngs <- make_KidProngs_from_FSL_and_PMMFS(ret$FSL, ret$PMMFS)
  
  # this is just an integer vector that will be a stack we push empty entries onto and off of
  ret$Pile <- integer(0)  # length(ret$IFS))
  
  # same with this
  ret$MatPile <- integer(0)
  
  # these are just the genotype frequencies
  ret$Gfreqs <- gfreqs_from_afreqs(Vars$afreqs) 
  
  # these are the unrelated pair genotype frequencies (3 x 3 x L)
  ret$UPG <- unrelated_pair_gfreqs(ret$Gfreqs)
  
  # these are the transition probs (3 x 3 x 3)
  ret$TP <- trans_probs()
  
  # these are the likelihoods of all the genotypes of the indivs in the sample
  ret$IndLiks <- Vars$geno_liks
  
  # now simulate to see what the distribution of logls is
  message("Staring LogL simulation at ", date(), "    This will take a few seconds")
  logls <- simulate_sib_pair_logls(ret$Gfreqs, 0.005, 100000)
  q1000 <- quantile(logls$LogL_ratio[logls$Relat=="Full_Sib_Pairs"], probs=.001)
  message("Done with LogL simulation at ", date())
  message("99.9% of true full sibs simulated to have logL > than ", format(q1000, digits=6))

  
  
  # this makes a list of vectors that have the base-0 indices of the possible
  # full sibs of each individual.  AFSL = Acceptable Full Sibling List
  ret$AFSL <- high_logl_pairs(FSP = as.vector(full_sibling_pair_gfreqs(ret$Gfreqs, mu)), 
                           UPF = as.vector(unrelated_pair_gfreqs(ret$Gfreqs)), 
                           G   = Vars$snp_genos$mat, 
                           loglV = q1000)
  
  message("\nDone with all mcmc variable initialization at ", date())
  
  ret  # return that big list
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