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


//' compute the genotype probs of a pair of full siblings with genotyping error
//' 
//' This could be done entirely in R, but it is hard to think it through and it seemed
//' it would be hard to maintain.  Super easy in C.
//' 
//' @param  L the number of loci
//' @param G the number of genotypic states
//' @param PP parent pair probs as returned by unrelated_pair_gfreqs() for example
//' @param TP transmision probs as returned by trans_probs()
//' @param GE prob of observed genotypes given true genotypes, as returned by lik_array_from_simple_geno_err() for example
//' @export
// [[Rcpp::export]]
NumericVector C_full_sibling_pair_gfreqs(int L, int G, NumericVector PP, NumericVector TP, NumericVector GE) {
  int l, o1, o2, a, b, t1, t2;
  NumericVector ret(G * G * L);
  double acc;  // for accumulating sums
  
  #define pp(a, b, l) PP[a + (b * G) + (l * G * G)]
  #define tp(a, b, c) TP[a + (G * b) + (G * G * c)]
  #define ge(t, o, l) GE[t + (G * o) + (G * G * l)]
  #define _ret(o1, o2, l) ret[o1 + (G * o2) + (G * G * l)]
  
  for(l=0; l<L; l++) {  // over loci
    for(o1=0; o1<G; o1++) {  // over sibling 1's observed genotype
      for(o2=0; o2<G; o2++) { // over sibling 2's observed genotype
        acc = 0.0;  // initialize to accumulate a sum
        for(a=0; a<G; a++) {  // over first parent geno
          for(b=0; b<G; b++) { // over second parent geno
            for(t1=0; t1<G; t1++) { // over sibling 1's true genotype
              for(t2=0; t2<G; t2++) { // over sibling 2's true genotype
                acc += pp(a, b, l) * tp(a, b, t1) * tp(a, b, t2) * ge(t1, o1, l) * ge(t1, o1, l) * ge(t2, o2, l);
              }
            }
          }
        }
        _ret(o1, o2, l) = acc;  // assign it here to the output
      }
    }
  }
  return(ret);
  
  #undef pp
  #undef tp
  #undef ge
  #undef _ret
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


//' Update one row of a marriage likelihood matrix IN PLACE!
//' 
//' The intended use of this is to update the marriage likelihoods after an individual has
//' been moved from one sibship to another.  So, for example, if you moved an individual 
//' from sibship 40 to sibship 291 (as subscripted in R) then you would want to run this
//' with \code{bz_idx} equal to \code{c(39, 290)}.  
//' @param S a list of vectors that give the indices (base 0) of the individuals in the 
//' full sibling groups.  
//' @param PK per-kid marriage likelihoods.  This must be of class \code{\link{marriage_geno_lik_array}},
//' which is just a matrix underneath with G x G x L rows and N columns.
//' @param ML the marriage likelihoods matrix to be modified
//' @param bz_idx An integer vector holding the BASE-0 indices of the components of S that 
//' will be accessed and used to update ML.
//' 
//' @return This doesn't return anything.  It modifies ML in place via call be reference.  Our
//' goal here is to make updates without copying a lot of memory.
//' @export
// [[Rcpp::export]]
void update_marriage_likelihoods_in_place(List S, NumericMatrix PK, NumericMatrix ML, IntegerVector bz_idx) {
  int yl;
  IntegerVector y;
  List tmp;
  int marriage;
  
  for(IntegerVector::iterator i = bz_idx.begin(); i != bz_idx.end(); ++i) {
    tmp = S[*i];
    y = as<IntegerVector>(tmp["Indivs"]);
    marriage = as<int>(tmp["LMMI_Idx"]);
    yl = y.length();
    ML( _, marriage) = PK( _, y[0]);   // this initializes it to accumulate a product
    for(int yi=1; yi<yl; yi++) {
      ML(_, marriage) =  ML(_, marriage) * PK(_, y[yi]);
    } 
  }
}



//' Calculate/update one row of a marriage Posterior matrix IN PLACE!
//' 
//' The intended use of this is to update the marriage posteriors after an individual has
//' been moved from one sibship to another.  So, for example, if you moved an individual 
//' from sibship 40 to sibship 291 (as subscripted in R) then you would want to run this
//' with \code{bz_idx} equal to \code{c(39, 290)}, after you had updated the LMMFS.  
//' @param S a list of lists each with two components.  The first is LMMI_Idx which is the base-0
//' index of the Marriage that the component is referring to, and the second is Indivs which give
//' the indices (base 0) of the individuals in the 
//' full sibling groups.  
//' @param PK per-kid marriage likelihoods.  This must be of class \code{\link{marriage_geno_lik_array}},
//' which is just a matrix underneath with G x G x L rows and N columns.
//' @param ML the marriage likelihoods matrix
//' @param MP the marriage posteriors matrix that is to be updated
//' @param Pri  the join prior probbabilities of the parent pair
//' @param NGS Number of genotypic states.  For pairs of parents, for example, this will be 9. 
//' @param bz_idx An integer vector holding the BASE-0 indices of the components of S that 
//' will be accessed and used to update ML.
//' 
//' @return This doesn't return anything.  It modifies ML in place via call be reference.  Our
//' goal here is to make updates without copying a lot of memory.
//' @export
// [[Rcpp::export]]
void update_marriage_posteriors_in_place(List S, NumericMatrix ML, NumericMatrix MP, NumericMatrix Pri, int NGS, IntegerVector bz_idx) {
  int yl;
  IntegerVector y;
  List tmp;
  int nL = ML.nrow() / NGS;  // this should be the number of loci
  int j,k,rr;
  double sum;
  int marriage;
  
  for(IntegerVector::iterator i = bz_idx.begin(); i != bz_idx.end(); ++i) {
    tmp = S[*i];
    marriage = as<int>(tmp["LMMI_Idx"]);
    MP( _, marriage) = ML( _, marriage) * Pri;   // multiply it by the prior
    
    // now we cycle over the loci and normalize the sums in each
    for(j=0; j<nL; j++) {
      sum = 0.0;
      for(k=0; k<NGS; k++) {
        rr = j * NGS + k;
        sum += MP(rr, marriage);
      }
      for(k=0; k<NGS; k++) {
        rr = j * NGS + k;
        MP(rr, marriage) /= sum;
      }
    }
  }
}




//' update marriage node kid prongs IN PLACE
//' 
//' This uses the information in the PMMFS matrix (the 
//' "Posterior Matrix of Marriages given Full Siblings") to calculate
//' the posterior predictive distribution for the next full sibling
//' in the full sibling group which is referenced through S and bz_idx
//'  
//' @param S a list of lists each with two components.  The first is LMMI_Idx which is the base-0
//' index of the Marriage that the component is referring to, and the second is Indivs which give
//' the indices (base 0) of the individuals in the 
//' full sibling groups.  
//' @param MP the marriage posteriors matrix to be used in the udpating. Should be a matrix
//' that has NGS_P * L rows, where L is the number of loci.
//' @param KP Kid prongs.  This is what gets updated in place. Should be a matrix that has
//' NGS_K * L rows, where L is the number of loci.
//' @param NGS_P Number of genotypic states in a parent pair.  This will typically be 9.
//' @param NGS_K Number of genotypic states in a kid.  This will typically be 3
//' @param TP transmision probs as returned by trans_probs()
//' @param bz_idx An integer vector holding the BASE-0 indices of the components of S that 
//' will be accessed and used to update ML.
//' 
//' @return This doesn't return anything.  It modifies Prongs in place via call by reference.  Our
//' goal here is to make updates without copying a lot of memory.
// [[Rcpp::export]]
void update_marriage_node_kid_prongs_in_place(List S, NumericMatrix MP, NumericMatrix KP, 
                                         int NGS_P, int NGS_K, NumericVector TP, IntegerVector bz_idx) {
  int yl;
  IntegerVector y;
  List tmp;
  int nL = KP.nrow() / NGS_K;  // this should be the number of loci
  int j,k, p;
  double sum;
  int marriage;
  
  #define KP_ KP(k + NGS_K * j, marriage)
  #define MP_ MP(p + NGS_P * j, marriage)
  #define TP_ TP(p + NGS_P * k)
  
  for(IntegerVector::iterator i = bz_idx.begin(); i != bz_idx.end(); ++i) {
    tmp = S[*i];
    marriage = as<int>(tmp["LMMI_Idx"]);
    for(j=0; j<nL; j++) {  // cycle over loci
      sum = 0.0;
      for(k=0; k<NGS_K; k++) {  // cycle over the 3 kid genotypes at each locus
        KP_ = 0.0;     // initialize to accumulate a su
        for(p=0; p<NGS_P; p++) { // cycle over the 9 parent genotypes at each locus
          KP_ += MP_ * TP_;
        }
        sum += KP_;
      }
      for(k=0; k<NGS_P; k++) {  // cycle over the kid genotypes one more time to normalize them to sum to one
        KP_ /= sum;
      }
    }
  }
  #undef KP_
  #undef MP_
  #undef TP_
}



// [[Rcpp::export]]
NumericVector kid_prongs_times_ind_likelihoods(List S, NumericVector IndGenoLik, NumericMatrix KidProngs) {
  List tmp;
  NumericVector ret(S.length());
  int fsp;
  double prod = 1.0;
  int j,k;
  int nL = IndGenoLik.length() / 3; 
  double sum;
  
  for(int i=0; i<S.length(); i++) {  // cycle over elements of S
    tmp = S[i];
    int fsp = as<int>(tmp["LMMI_Idx"]);
    prod = 1.0;  // prepare to accumulate a product
    
    for(j=0; j<nL; j++) {
      sum = 0.0;
      for(k=0; k<3; k++) {
        sum+= KidProngs(k + 3 * j, fsp) * IndGenoLik[k + 3 * j];
      }
      prod *= sum;
    }
    ret[i] = prod;
  }
  return ret;
}


// [[Rcpp::export]]
List gibbs_update_one_indiv_in_place( List FSL,
                                      IntegerVector IFS,        
                                      NumericMatrix LMMI,       
                                      NumericMatrix LMMFS,      
                                      NumericMatrix PMMFS,      
                                      NumericMatrix KidProngs,  
                                      IntegerVector Pile,       
                                      NumericMatrix Gfreqs,     
                                      NumericMatrix UPG,        
                                      NumericVector TP,         
                                      NumericMatrix IndLiks,    
                                      int Ind                   
                                      ) {
    int fs = IFS[Ind];
    List tmp = FSL[fs];
    int fsp = tmp["LMMI_Idx"];  // this is the column that Ind's current full sib group occuppies in LMMFS
                                // and PMMFS, and KidProngs
    IntegerVector current_sibs = tmp["Indivs"];  // everyone in Ind's current sibgroup (including himself)
    int n = current_sibs.length();
    List ret(5);
                                
    // now make some copies of what it looked like before the update (easier to go back to this)                            
    NumericVector LMMFS_orig_col = LMMFS( _, fsp);
    NumericVector PMMFS_orig_col = PMMFS( _, fsp);
    NumericVector KidProngs_orig_col = KidProngs( _, fsp);
    
    if(n==1) {  // if Ind is the only one in his current sibgroup, then use that sibgroup to see if he 
                // wants to be a solo sib-group still.
      KidProngs( _, fsp) = Gfreqs; 
    }
    else if(n>1) { // if there are others in his sibgroup, divide out his likelihood contribution to the 
                   // the LMMFS and recompute PMMFS and KidProngs
       LMMFS( _, fsp) = LMMFS( _, fsp) / LMMI( _, Ind);
       update_marriage_posteriors_in_place(FSL, LMMFS, PMMFS, UPG, 9, fs);
       update_marriage_node_kid_prongs_in_place(FSL, PMMFS, LMMI, 9, 3, TP, fs);
    }
    else  {
      ; // need to throw an error here.
    }
    
    
    // now, we can zoom over all the KidProngs and compute the full conditional likelihood of 
    // Ind belonging to each. 
    NumericVector FCLs = kid_prongs_times_ind_likelihoods(FSL, IndLiks( _, Ind), KidProngs);
    
    ret[0] = FCLs;
    
    return(ret);
}

