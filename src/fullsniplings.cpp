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
  NumericVector ol_dims = offspring_likelihoods.attr("gqa_dim");
  List ol_dn = offspring_likelihoods.attr("gqa_dimnames");
  int G = ol_dims[0],  // G is the # of genotypes.  It should always be 3, but i'll let it vary
      L = ol_dims[1],  // Number of loci
      N = ol_dims[2];  // Number of individuals 
  NumericVector ret(G * G * L * N);  // initialized to 0's with length G*G*L*N
  NumericVector rd = NumericVector::create(G, G, L, N);  // dimensions of output. Though we will squash it into matrix in the next line
  ret.attr("dim") = NumericVector::create(G * G * L, N);
  
  // set all the attributes to make it a marriage_geno_lik_array object
  ret.attr("gqa_dim") = rd;
  ret.attr("gqa_dimnames") = List::create( _["Parent1_Geno"] = ol_dn[0],
                                            _["Parent2_Geno"] = ol_dn[0],
                                            _["Loci"]         = ol_dn[1],
                                            _["Offspring"]    = ol_dn[2]  
                                            );
  ret.attr("class") = CharacterVector::create("marriage_geno_lik_array", "geno_qty_array");
  ret.attr("gqa_description") = CharacterVector::create("parent pair likelihoods given each, single offspring");
  
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


//' compute the marriage node likelihoods given offspring specified in a full-sibling list
//' 
//' @param S a list of vectors that give the indices (base 0) of the individuals in the 
//' full sibling groups.  
//' @param PK per-kid marriage likelihoods.  This must be of class \code{\link{marriage_geno_lik_array}},
//' which is just a matrix underneath with G x G x L rows and N columns.
//' @export
// [[Rcpp::export]]
NumericMatrix multi_kid_marriage_likelihoods(List S, NumericMatrix PK) {
  int M = S.length();  // number of marriages we are computing for here
    IntegerVector PD = PK.attr("gqa_dim");
    List PKdn = PK.attr("gqa_dimnames");
    int G1 = PD[0];
    int G2 = PD[1];
    int L  = PD[2];
    int m=0, yl;
    IntegerVector y;
    NumericMatrix ret(G1 * G2 * L, M); // allocate memory to be returned
    
    
   
    // cycle over the elements of M
  for(List::iterator it=S.begin(); it != S.end(); ++it)  {
    y = as<IntegerVector>(*it);
     yl = y.length();
     ret( _, m) = PK( _, y[0]);   // this initializes it to accumulate a product
     for(int yi=1; yi<yl; yi++) {
       ret(_, m) =  ret(_, m) * PK(_, y[yi]);
     } 
     m++;
    }
   
    // set attributes and make it of class marriage_geno_lik_array
    ret.attr("gqa_dim") = NumericVector::create(G1, G2, L, M);
    ret.attr("class") = CharacterVector::create("marriage_geno_lik_array", "geno_qty_array");
    ret.attr("gqa_description") = CharacterVector::create("parent pair likelihoods given all the offspring currently allocated to them");
    ret.attr("gqa_dimnames") = List::create(_["Parent1_Geno"] = PKdn[0],
                                            _["Parent2_Geno"]   = PKdn[1],
                                            _["Loci"]           = PKdn[2],
                                            _["Marriages"]      = S.attr("names")
                                            );
   return(ret);
}

