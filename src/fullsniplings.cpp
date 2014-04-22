#include <Rcpp.h>
using namespace Rcpp;




//' Compute the likelihood of the genotypes of a marriage given the genotype of \emph{just one} of its offspring
//' 
//' It is advantageous to precompute and store some quantities for pedigree analysis with SNPs.
//' One of those quantities is the likelihood of the genotypes of the two parents that produced a child,
//' conditional only on that one child's genotype. This is convenient because the likelihood of the two parents' 
//' genotypes given \emph{all} of their offspring will be the product of the single-offspring
//' likelihoods.  When genotyping error is present, that has to be thrown into the mix as well.
//' However, here we assume that the individual offspring genotype likelihoods have already been
//' computed given their observed genotypes.
//' @param offspring_likelihoods This is a 3 x L x N array of likelihoods of the individual's
//' true genotypes given their observed genotypes.  This sort of array is returned by the function
//' \code{\link{get_indiv_geno_lik}}.  
//' @param transmission_probs a 3 x 3 x 3 matrix of probs of kid genotypes given parent genotypes.
//' Such an array is returned by \code{\link{trans_probs}}.  See the documentation thereof for a
//' description.
//' @return This returns a 3 x 3 x L x N array of the parent genotypes (9 states) at each of the L loci
//' for each of the N offspring.  Although these are ostensible likelihoods of "marriages" (think of them
//' as marriage nodes), it will be beneficial to regard these quantities of properties of the individual
//' offspring moreso than as properties of the marriage nodes.
//' @export
// [[Rcpp::export]]
NumericVector per_kid_marriage_likelihoods(NumericVector offspring_likelihoods, NumericVector transmission_probs) {
  NumericVector ol_dims = offspring_likelihoods.attr("dim");
  int G = ol_dims[0],  // G is the # of genotypes.  It should always be 3, but i'll let it vary
      L = ol_dims[1],  // Number of loci
      N = ol_dims[2];  // Number of individuals 
  NumericVector ret(G * G * L * N);  // initialized to 0's with lenth G*G*L*N
  NumericVector rd = NumericVector::create(G, G, L, N);  // dimensions of output.
  ret.attr("dim") = rd;
  
// The following error catching just crashes my machine.  Leave it out for now.
// I don't know what is up with Rcpp on that.
//  if(ol_dims.size() != 3) {
//    throw Rcpp::exception("Unexpected condition occurred");
//  }
//  if(ol_dims[0] != 3) stop("offspring likelihoods first dimension must be of length 3");
 

  // for indexing into the offspring_likelihoods array
  #define IDX(genotype, locus, indiv)  genotype + (locus * G) + (indiv * G * L) 
  
  // for indexing into the ret array
  #define IDX2(par1, par2, locus, indiv) par1 + (par2 * rd[0]) + (locus * rd[0] * rd[1]) + (indiv * rd[0] * rd[1] * rd[2])
  
  // for indexing into the transmission_probs array
  #define IDX3(par1, par2, kid)    par1 + (G * par2) + (G * G * kid)
  
  for(int ind=0; ind<N; ind++) {
    for(int loc=0; loc<L; loc++) {
      for(int pa=0; pa<G; pa++) {
        for(int ma=0; ma<G; ma++) {
          // at this level we will do calculations for the pa, ma, genotype combo at locus loc for the parents of kid ind
          for(int kg=0; kg<G; kg++) {  // kg is "kid geno"
            ret[IDX2(pa, ma, loc, ind)] += transmission_probs[IDX3(pa, ma, kg)] * offspring_likelihoods[IDX(kg, loc, ind)];
          }
        }
      }
    }
  }


  return(ret);
} 



//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x, NumericVector y) {
   return x * y * 2;
}
