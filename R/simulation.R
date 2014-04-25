

#' simulate pairs of full sibs and pairs of unrelateds and compute loglikelihoods
#' 
#' @param gf a genotype frequencies array
#' @param Reps the number of replicates to draw
#' @param Genotyping error rate (gets recycled over loci as necessary)
#' @export
simulate_sib_pair_logls <- function(gf, mu=.005, Reps=10000) {
  SP_mat <- full_sibling_pair_gfreqs(gf, mu)
  UR_mat <- unrelated_pair_gfreqs(gf)
  if(ncol(SP_mat) != ncol(UR_mat)) stop("Column number mismatch between SP_mat and UR_mat within function.")

  SS <-lapply(1:ncol(SP_mat), function(x){
      list(
        sp_samp = sample(1:nrow(SP_mat), Reps, replace=T, prob = SP_mat[, x]),
        ur_samp = sample(1:nrow(SP_mat), Reps, replace=T, prob = UR_mat[, x]),
        sp_logprobs = log(SP_mat[, x]),
        ur_logprobs = log(UR_mat[, x])
      )
    })
  
  # now, a simple vanilla MCMC would be like this
  ret <- stack(list(
    Full_Sib_Pairs =   unname(rowSums(sapply(SS, function(x) x$sp_logprobs[x$sp_samp] - x$ur_logprobs[x$sp_samp]))),
    Unrelated_Pairs = unname(rowSums(sapply(SS, function(x) x$sp_logprobs[x$ur_samp] - x$ur_logprobs[x$ur_samp])))
  ))
  
  names(ret) <- c("LogL_ratio", "Relat")
  ret
}

