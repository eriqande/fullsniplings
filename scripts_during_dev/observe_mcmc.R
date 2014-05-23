## Script to test and observe gibbs update steps

##  This started out in full-sibling-chinook-mcmc.Rmd but then I put it here.

library(fullsniplings)
library(pryr)


# this is a function to summarize what happened/changed while running the function:
# x is the list that gibbs_update spits out
summarize_results <- function(x) {
  
  # changes in IFS
  IFS_changes = cbind(1:length(x$IFS_start), x$IFS_start, x$IFS_end)[x$IFS_start != x$IFS_end, ]
  if(length(IFS_changes)==3) names(IFS_changes) <- c("Ind", "Start", "End")
  
  # here we identify the changes in FSL
  differs <- which(sapply(1:length(x$FSL_start), function(y) !identical(x$FSL_start[[y]], x$FSL_end[[y]])))
  
  FSL_changes <- lapply(differs, function(j) paste("FS ", j, "  ", "LMMI_Idx:  ", x$FSL_start[[j]]$LMMI_Idx, " --> ", x$FSL_end[[j]]$LMMI_Idx,
                                    "         Indivs:  ", paste(x$FSL_start[[j]]$Indivs, collapse=" "), "  -->  ", paste(x$FSL_end[[j]]$Indivs, collapse=" "),
                                    sep=""))
  
  
  list(
    IFS_changes = IFS_changes,
    FSL_changes = FSL_changes
    )
  
}




# here is a function that susses out whether the high likelihood sibships have full siblings in them
# or not
# x is the list returned by gibbs_update
# sgl is the sg.list
comp_liks_to_num_true_sibs <- function(x, sgl) {
  
  # here are the LMMI_Idx's of all those full sibling groups too
  lmmi_idx <- sapply(x$FSL_start[x$AFS+1], function(y) y$LMMI_Idx)
  
  df <- data.frame(AFS = x$AFS, LMMI_Idx = lmmi_idx,  FCLs = x$FCLs)  # ultimately we will return this thing
  
  # these are the individuals in these full-sibling groups
  fsl <- lapply(x$FSL_start[df$AFS+1], function(y) y$Indivs)
  
  
  
  df$NumSibs <- sapply(fsl, function(z) length(intersect(z, sgl[[x$Ind+1]])))
  
  df$SibIndices <- sapply(fsl, function(z) paste(intersect(z, sgl[[x$Ind+1]]), collapse=",") )
  
  df[order(df$FCLs, decreasing = TRUE), ]
  
}


# this function tables the sibship sizes
table_sibsizes <- function(S) {
  tab <- table(sapply(S, function(x) length(x$Indivs)))
  tot <- sum(as.numeric(names(tab)) * tab)
  list(Table=tab, Total=tot)
}

# this function orders the FSL by size and orders the indivs within each FSL for easy viewing
ordered_fsl <- function(S) {
  inds <- lapply(S, function(x) sort(x$Indivs))
  ord <- order(sapply(inds, length), decreasing=T)
  inds[ord]
}


# this function returns a vector of strings naming the current sibships, like 
# 10-34-67-98, which can be used to hash them
hashable_sibships <- function(S) {
  osibs <- ordered_fsl(S)   # ordered sibgroups
  osibs <- osibs[sapply(osibs, length)>0]  # drop the empty ones
  sapply(osibs, function(x) paste(x, collapse="-"))
}


# also, before starting, it will be nice to make a list for every individual with the base-0 subscripts of the other
# individuals in his full sibling group:
ped <- fs_dev_test_data$chinook_full_sibs_pedigree
ordered.kids <- rownames(fs_dev_test_data$chinook_full_sibs_genos)
ped$kididx <- as.integer(factor(ped$Kid, levels = ordered.kids)) - 1
kid.split <- split(ped$kididx, paste(ped$Pa, "---", ped$Ma, sep=""))

# here we get a list of all the sibling groups
sibgroups <- unname(kid.split[order(sapply(kid.split, length), decreasing=T)])

# and here we can store them as strings:
true_sibgroups_as_strings <- sapply(sibgroups, function(x) paste(x, collapse="-"))


# now get a list of the siblings of an individual
sg.list <- vector("list", length = nrow(ped))
for(L in kid.split) {
  for(ind in L) {
    sg.list[[ind+1]] <- setdiff(L, ind)
  }
}
# If Ind is the base-0 index of an individual then
# sg.list[[Ind+1]] is a vector of the base-0 indices of all of its full siblings




# initialize the chain to everyone in their own sibgroup.
#chinook_chain <- full_sib_mcmc_initialize(fs_dev_test_data$chinook_full_sibs_genos, mu = 0.005)
#save(chinook_chain, file="chinook_chain.rda")
load("chinook_chain.rda")  # i created the rda with the above two and then just do this to save time

# make a hard copy of that
cc <- chinook_chain

cc$PMMFS <- unclass(cc$PMMFS)

# let's see how many of the total number of siblings do *not* appear in the acceptable full
# sibling list
afs_check <- t(sapply(seq_along(sg.list), function(x) c(length(sg.list[[x]]), length(sg.list[[x]] %in% cc$AFSL[[x]]) )))
colnames(afs_check) <- c("NumSibs", "NumSibsInAFS")
afs_check[ (afs_check[,1] - afs_check[,2]) != 0, ]  # Far out.  We aren't missing any sibs in the AFSL.



# here is the basic way we run this:
if(FALSE) {
set.seed(5)
i <- 248
grab <- do.call(what = gibbs_update_one_indiv_in_place, args = (c(cc, Ind=i)))
cc$Pile <- grab$Pile  # currently these have to be copied back
cc$MatPile <- grab$MatPile


summarize_results(grab)
grab$solo_lik
head(comp_liks_to_num_true_sibs(grab, sg.list), n=20)
}




# here we can run it multiple times on a specific group of individuals:
set.seed(5)
#for(j in 1:40) {for(i in sibgroups[[1]]) {
visited_sibgroups <- vector(mode = "integer")
burn_in <- 50
num_sweeps <- 100
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
  }
}



# this is parallel to visited_sibgroups
vis_sib_lengths <- sapply(strsplit(names(visited_sibgroups), "-"), length)
names(vis_sib_lengths) <- names(visited_sibgroups)

# extract all the true sibgroups from there and see how many times they were visited (or NA if they were not visited)
visited_sibgroups[true_sibgroups_as_strings]

# here we see how many of our 100% posterior sibships were correct
names(visited_sibgroups)[visited_sibgroups==100] %in% true_sibgroups_as_strings

# so, what would be really nice is a function that summarizes this and gives us the sizes
# of the correct and the sizes of the incorrect sibships.  Hmmm...we should be able
# to put all this into "long format" with a few columns.  Say:
#   SibString(as rownames)   SibSize     Posterior     Correct    
run_results <- data.frame(Posterior = visited_sibgroups)
run_results$SibSize <- vis_sib_lengths
run_results$Correct <- names(visited_sibgroups) %in% true_sibgroups_as_strings
# and then order it by posterior and then size
run_results <- run_results[order(run_results$Posterior, run_results$SibSize, decreasing=T), ]

# now get the cumulative number of individuals placed into the correct sibships (and the wrong one)
run_results$CusumCorrect <- cumsum(run_results$SibSize * run_results$Correct)
run_results$CusumInCorrect <- cumsum(run_results$SibSize * run_results$Correct==FALSE)

# Note! I really need to to something different for the above.  Maybe a greedy approach where I take every sibship
# as I scan down the list, and store which individuals I have seen, and then don't accept any sibgroups containing
# individuals I have already inferred to be in higher-posterior sibgroups.  

# and I should also compute partition distances.

# It is looking good though!


