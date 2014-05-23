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
  df <- data.frame(AFS = x$AFS, FCLs = x$FCLs)  # ultimately we will return this thing
  
  # these are the individuals in these full-sibling groups
  fsl <- lapply(x$FSL_start[df$AFS+1], function(y) y$Indivs)
  
  df$NumSibs <- sapply(fsl, function(z) length(intersect(z, sgl[[x$Ind+1]])))
  
  df$SibIndices <- sapply(fsl, function(z) paste(intersect(z, sgl[[x$Ind+1]]), collapse=",") )
  
  df[order(df$FCLs, decreasing = TRUE), ]
  
}
 

# also, before starting, it will be nice to make a list for every individual with the base-0 subscripts of the other
# individuals in his full sibling group:
ped <- fs_dev_test_data$chinook_full_sibs_pedigree
ordered.kids <- rownames(fs_dev_test_data$chinook_full_sibs_genos)
ped$kididx <- as.integer(factor(ped$Kid, levels = ordered.kids)) - 1
kid.split <- split(ped$kididx, paste(ped$Pa, "---", ped$Ma, sep=""))
sg.list <- vector("list", length = nrow(ped))
for(L in kid.split) {
  for(ind in L) {
    sg.list[[ind+1]] <- setdiff(L, ind)
  }
}

# cool, at this point, if Ind is the base-0 index of an individual then
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
set.seed(5)
i <- 955
grab <- do.call(what = gibbs_update_one_indiv_in_place, args = (c(cc, Ind=i)))
cc$Pile <- grab$Pile  # currently these have to be copied back
cc$MatPile <- grab$MatPile

summarize_results(grab)
grab$solo_lik
head(comp_liks_to_num_true_sibs(grab, sg.list), n=20)

