


#' Run the MCMC for full sibling inference
#' 
#' @inheritParams full_sib_mcmc_initialize
#' @param burn_in Number of sweeps to discard as burn-in.  (In one sweep, every
#' individual in the sample has been updated once by Gibbs sampling.)
#' @param num_sweeps Number of sweeps after burn in to use as a sample for 
#' MCMC.
#' @export
#' @examples
#' # this is a much shorter run than one should do:
#' run_mcmc(fs_dev_test_data$chinook_full_sibs_genos, burn_in = 3, num_sweeps = 5)
run_mcmc <- function (genos, mu = 0.005, pair_prob_cutoff = 0.001, burn_in, num_sweeps) {
  # initialize a chain and call it cc
  cc <- full_sib_mcmc_initialize(genos, mu, pair_prob_cutoff)
  
  # set up a vector to store a hash
  visited_sibgroups <- vector(mode = "integer")
  for(j in 1:(burn_in + num_sweeps)) {
    for(i in 0:(length(cc$IFS)-1)) {
      grab <- do.call(what = gibbs_update_one_indiv_in_place, args = (c(cc, Ind=i)))
      cc$Pile <- grab$Pile  # currently these have to be copied back
      cc$MatPile <- grab$MatPile
    }
    print(table_sibsizes(S = cc$FSL)$Table)
    if(j>burn_in) {
      sibstrs <- hashable_sibships(cc$FSL)
      visited_sibgroups[sibstrs][is.na(visited_sibgroups[sibstrs])] <- 0  # set them to 0 if they were not visited previously
      visited_sibgroups[sibstrs] <- visited_sibgroups[sibstrs] + 1
      print(paste("Completed sample collection sweep", j-burn_in, "of", num_sweeps, "at", date()))
    }
    else print(paste("Completed burn-in sweep", j, "of", burn_in, "at", date()))
  }
  visited_sibgroups
  
  GP <- greedy_partition(visited_sibgroups)
  Partition <- as.data.frame(t(sapply(GP$VSL, function(x) c(x$Visits, x$NumSibs))))
  names(Partition) <- c("Visits", "NumSibs")
  SibGroupsByName <- lapply(strsplit(rownames(Partition), "-"), function(x) rownames(genos)[as.numeric(x) + 1])
  list(Partition = Partition, GP_result = GP, VS = sort_visited_sibgroups(visited_sibgroups), 
       SibGroupsByName = SibGroupsByName)
}


#' sorts the visited sibgroups on the basis of posterior and size
#' 
#' Returns a list that has for each sibgroup a list which holds: 
#' a vector of \code{Indivs} the number of visits, \code{Visits},
#' and the size of the sibgroup, \code{NumSibs}
#' @param VS A vector of counts of the number of times each sibship was visited.  The sibships are named
#' according to who is in them via the hashable_sibships function.
sort_visited_sibgroups <- function(VS) {
  # make a list that includes components Visits (for number of times visited) and also 
  # a vector of the sibling groups
  VSL <- lapply(names(VS), function(x) list(Visits = VS[x], Indivs = as.integer(strsplit(x, "-")[[1]])))
  names(VSL) <- names(VS)
  VSL <- lapply(VSL, function(x) c(x, NumSibs = length(x$Indivs)))  # add the number of sibs to the list
  
  VSL <- VSL[order(sapply(VSL, function(x) x$Visits), sapply(VSL, function(x) x$NumSibs), decreasing=c(T,F))]
  VSL
}

#' Uses the estimated posteriors for each sibgroup to come up with a partition into sibgroups
#'
#' Used internally.  This is a quick and dirty greedy thing. It basically orders sibships 
#' by posterior and then by size, and then it goes through and calls them sibships starting with
#' the highest posterior ones.  If it reaches one that includes someone that is already in a sibship
#' it discards it and continues.  At the end, any one not assigned to a sibship is made a singleton.
#' @inheritParams run_mcmc
#' @inheritParams sort_visited_sibgroups
greedy_partition <- function(VS) {
  
  VSL <- sort_visited_sibgroups(VS)
  Dumped <- VSL # for storing info about who gets dumped and why
  
  Everyone <- sort(unique(unlist(lapply(VSL, function(x) x$Indivs))))
  # now, we cycle over the elements of VSL and toss those elements that have individuals that 
  # were already seen
  got_em <- integer(0)
  for(i in seq_along(VSL)) {
    Intersecters <- intersect(VSL[[i]]$Indivs, got_em)
    if(length(Intersecters) == 0) {  # if we haven't seen any of these sibs before
      got_em <- c(got_em, VSL[[i]]$Indivs)  # add them to got_em and set their entry to a 0-length vector
      Dumped[[i]] <- integer(0)
    }
    else { # if we have seen them.  NA them out of VSL and say why in Dumped.
      VSL[[i]] <- integer(0)
      Dumped[[i]]$Intersecters <- Intersecters
    }
  }
  
  VSL <- VSL[sapply(VSL, function(x) length(x) > 0)]
  Dumped <- Dumped[sapply(Dumped, function(x) length(x) > 0)]
  got_em <- sort(got_em)
  missed_em <- setdiff(Everyone, got_em)
  tmp <- lapply(missed_em, function(x) list(Visits = 0, Indivs = x, NumSibs = 1)) # put these in here with zero visits
  names(tmp) <- missed_em
  VSL <- c(VSL, tmp)
  
  list(VSL = VSL, Dumped = Dumped)
  
}
