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
                acc += pp(a, b, l) * tp(a, b, t1) * tp(a, b, t2) * ge(t1, o1, l) * ge(t2, o2, l);
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
//' with \code{bz_idx} equal to the column of sibship 39, and then run it again with it equal 
//' to the column of  290, after you had updated the LMMFS
//' for each of those.
//' @param ML the marriage likelihoods matrix
//' @param MP the marriage posteriors matrix that is to be updated
//' @param Pri  the joint prior probbabilities of the parent pair
//' @param NGS Number of genotypic states.  For pairs of parents, for example, this will be 9. 
//' @param bz_idx An integer vector holding the BASE-0 index of the columns of MP that needs updating.  Note that this
//' is not the index of the sibhip, but rather the column that sibship occupies in MP.  That will be something like 
//' FSL[Indiv]$LMMI_Idx that 
//' 
//' @return This returns the value of the index of the last marriage that it updated.  It does this so I can check some things.  It modifies MP in place via call by reference.  Our
//' goal here is to make updates without copying a lot of memory.
//' @export
// [[Rcpp::export]]
int update_marriage_posteriors_in_place(NumericMatrix ML, NumericMatrix MP, NumericMatrix Pri, int NGS, IntegerVector bz_idx) {
  IntegerVector y;
  int nL = ML.nrow() / NGS;  // this should be the number of loci
  int j,p;
  double sum;
  int marriage;
  
  // in the #defines below, p is for the 9 parent-genotype states
  #define MP_ MP[(NGS * nL * marriage) + (NGS * j)  +  p]
  #define ML_ ML[(NGS * nL * marriage) + (NGS * j)  +  p]
  #define Pri_ Pri[(NGS * j) + p]
  
  
  for(IntegerVector::iterator i = bz_idx.begin(); i != bz_idx.end(); ++i) {
    marriage = *i;
    
    // create the MP values by multiplying ML by the Prior
    for(j=0; j<nL; j++) {
      for(p=0; p<NGS; p++) {
        MP_ = ML_ * Pri_;
      }
    }
    
    // now we cycle over the loci and normalize the sums in each
    for(j=0; j<nL; j++) {
      sum = 0.0;
      for(p=0; p<NGS; p++) {
        sum += MP_;
      }
      for(p=0; p<NGS; p++) {
        MP_ /= sum;
      }
    }
  }
  
  #undef MP_
  #undef ML_
  #undef Pri_
  
  return(marriage);
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



//' Computes the posterior predictive likelihoods that an individual with likelihood vector IndGenoLik belongs to each sibling group in AFS
//' 
//' Note that we have to pass in the whole full sibling list (FSL) so that we can get the LMMI_Idx for each sibling group.
//' @param IndGenoLik Likelihoods of the focal individual.  There will be 3*L elements in that.  Note that \code{AFS} is the acceptable
//' sibling list for this focal individual.
//' @export
// [[Rcpp::export]]
NumericVector kid_prongs_times_ind_likelihoods(List FSL, NumericVector IndGenoLik, NumericMatrix KidProngs, IntegerVector AFS) {
  List tmp;
  int n = AFS.length();  // the number of sibgroups in AFS
  NumericVector ret(n);  
  double prod = 1.0;
  int j,k;
  int nL = IndGenoLik.length() / 3; 
  double sum;
  
  for(int i=0; i<n; i++) {  // cycle over the acceptable sibling groups (elements in AFS)
    tmp = FSL[AFS[i]];  // for each sibgroup, we need to ge its LMMI_Idx
    int fsp = as<int>(tmp["LMMI_Idx"]);  // this is its column in the KidProngs matrix
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





//' return a vector of indices of sibships that include at least one hi-sibship-lod individual
//' 
//' Note that it is currently up to the user to ensure that no element of AFS exceeds the length-1
//' of IFS
//' 
//' @param IFS the individual full siblings group vector. A vector such that element i contains the 
//' base-0 index of the sibship to which the individual with index i belongs.
//' @param AFS the vector of individual indexes (base 0) that have hi sibship lod with the 
//' focal individual.  i.e. AFS = the Acceptable Full Siblings
//' @export
//' @examples
//' possible_sibgroups(rep(16,50), 0:49)
//' possible_sibgroups(rep(0:9, each=3), sample(0:29, 10))
// [[Rcpp::export]]
IntegerVector possible_sibgroups(IntegerVector IFS, IntegerVector AFS) {
  int n = AFS.length();
  IntegerVector sg(n);  // to store sibgroup indices
  
  for(int i=0; i<n; i++) {
    sg[i] = IFS[AFS[i]];
  }
  
  return(sort_unique(sg));  // return each sibgroup just once
}




//' compute the pseudo-prior for individual and individual pulled out of sibgroup IndG
//' 
//' @param FSL the full sibling list
//' @param IndG The full sibling group to which the focal individual belongs
//' @param AFS acceptable full sibling groups for the focal individual.
//' @return This returns a list of two components. $solo is the prob that the individual will form a singleton and 
//' $afs is the prior that the individual will join any of the existing sibgroups in afs, scaled to sum to one.
//' @export
// [[Rcpp::export]]
List pseudo_prior(List FSL, int IndG, IntegerVector AFS) {
  List ret(2);
  int fn = FSL.length();
  int n = AFS.length();
  IntegerVector counts(fn);
  int NumOnes = 0;
  int Tot = 0;
  
  
  for(int i=0; i<fn; i++) {
    List tmp = FSL[i];
    IntegerVector Indivs = tmp["Indivs"];
    counts[i] = Indivs.length();
    if(i == IndG) counts[i] -= 1;  // remove him from his own sibgroup
    if(counts[i] == 1) NumOnes++;
    Tot += counts[i];
  }
  ret[0] = counts;
  ret[1] = NumOnes;
  
  NumericVector res(n);
  double sum = 0.0;
  for(int i=0; i<n; i++) {
    res[i] = counts[AFS[i]];
    sum += res[i];
  }
  for(int i=0; i<n; i++) {
    res[i] /= sum;
  }
  
  double solo_prob = (NumOnes * 1.0)/Tot;
  if(solo_prob > .995) solo_prob = .995;  // a quick hack to deal with the starting situation when everyone is a singleton.
  
  
  return List::create(_["solo"] = solo_prob, _["afs"] = res);
}




//' compute the posterior of Inds genotype given the genotype freqs in the pop
//' 
//' @param Gfreqs  array of genotype frequencies
//' @param Liks  array of likleihoods (3 * L in length)
//' @export
// [[Rcpp::export]]
double geno_post_c(NumericVector Gfreqs, NumericVector Liks) {
  int n = Liks.length();
  int L = n/3;
  double sum;
  double prod = 1.0;
  
 
  for(int i=0; i<L; i++) {
    sum=0.0;
    for(int j=0; j<3; j++) {
      sum += Gfreqs[i*3 + j] * Liks[i*3 + j];
    }
    prod *= sum;
  }
  return prod;
}


//' does a gibbs update of the full sibling group of individual Ind
//' @param FSL list of integer vectors. Each component is a list of two componontes: \code{LMMI_Idx} = the index of the sibgroup
//' in LMMI (see below).  It also is the column that Ind's current full sib group occuppies in LMMFS
//' @param IFS individual full siblings group vector. A vector such that element i contains the 
//' base-0 index of the sibship to which the individual with index i belongs.
//' and PMMFS, and KidProngs \code{Indivs} = the the base-0 indices of individuals in the full sibship.
//' @param LMMI Likelihood matrix of marriages given individuals.
//' @param LMMFS Likelihood matrix of marriages given full sibships.
//' @param PMMFS Posterior matrix of marriages given full sibships.
//' @param KidProngs Posterior predictives for the next individuals to be sampled from a full sibship.
//' @param Pile  An integer vector that we will use as a stack to hold the indices of empty FSL list elements
//' that we can put newly formed sibships into
//' @param MatPile An integer vector that is parallel to Pile that holds the corresponding column of LMMFS that 
//' go along with the elements of Pile. 
//' can be renewed.
//' @param AFSL The acceptable full siblings list.  Component i is a vector of the base-0 indices of the individuals
//' that indiv i has high full-sibling LOD with.
//' @param Gfreqs  The genotype frequencies.
//' @param UPG Unrelated Pair Genotype Frequencies.
//' @param TP Transmission probabilities
//' @param IndLiks  Individual likelihoods
//' @param Ind the index of the individual to be updated.  
//' 
//' 
//' @export
// [[Rcpp::export]]
List gibbs_update_one_indiv_in_place( List FSL,
                                      IntegerVector IFS,        
                                      NumericMatrix LMMI,       
                                      NumericMatrix LMMFS,      
                                      NumericMatrix PMMFS,      
                                      NumericMatrix KidProngs,  
                                      std::vector<int> Pile, 
                                      std::vector<int> MatPile,
                                      List AFSL,
                                      NumericMatrix Gfreqs,     
                                      NumericMatrix UPG,        
                                      NumericVector TP,         
                                      NumericMatrix IndLiks,    
                                      int Ind                   
                                      ) {
  
  // right at the top, we create a list to return values so we can see how things are progressing while developing this
  List ret = List::create(
      _["Ind"] = Ind,
      _["fs"] = NA_INTEGER,
      _["fsp"] = NA_INTEGER,
      _["fsp_vec"] = IntegerVector::create(NA_INTEGER),
      _["IFS_start"] = IntegerVector::create(NA_INTEGER),
      _["IFS_end"] = IntegerVector::create(NA_INTEGER),
      _["Pile"] = IntegerVector::create(NA_INTEGER),
      _["MatPile"] = IntegerVector::create(NA_INTEGER),
      _["FSL_start"] = clone(FSL),
      _["FSL_end"] = List(1),
      _["AFS"] = IntegerVector::create(NA_INTEGER),
      _["LMMFS_orig"] = NumericVector::create(NA_REAL),
      _["LMMFS_test"] = NumericVector::create(NA_REAL),
      _["PMMFS_orig"] = NumericVector::create(NA_REAL),
      _["PMMFS_test"] = NumericVector::create(NA_REAL),
      _["KidProngs_orig"] = NumericVector::create(NA_REAL),
      _["KidProngs_test"] = NumericVector::create(NA_REAL),
      _["FCLs"] = NA_REAL,
      _["solo_lik"] = NA_REAL,
      _["MP_Update_Idx1"] = NA_INTEGER
    );
    
  
  RNGScope Scope;  // Initialize random number generator.
  
  
ret["IFS_start"] = clone(IFS);
    int fs = IFS[Ind];   // index of the full sib group that Ind currently belongs to
    IntegerVector fs_vec = IntegerVector::create(fs);
ret["fs"] = fs;
    List tmp = FSL[fs];  
    int fsp = as<int>(tmp["LMMI_Idx"]);  // this is the column that Ind's current full sib group occuppies in LMMFS
                                // and PMMFS, and KidProngs
    IntegerVector fsp_vec = IntegerVector::create(fsp);  // this is messed up.  I need to do this to pass it to the functions that update MP's and KidProngs.
ret["fsp"] = fsp;
ret["fsp_vec"] = fsp_vec;
    IntegerVector current_sibs = tmp["Indivs"];  // everyone in Ind's current sibgroup (including himself)
    int n = current_sibs.length();  // how many individuals in the full sibgroup that Ind belongs to
    int New_fs;  // for storing the index of the new sibship and individual will be put into
    int New_fsp; // for storing the column of the LMMFS matrix that a new sibship will use.  
    IntegerVector New_fs_vec;
    IntegerVector New_fsp_vec;
    // Get the sibships that we would consider adding this individual into:
    // ASG = Acceptable Sibling Groups
    IntegerVector AFS = possible_sibgroups(IFS, AFSL[Ind]);  // note that this will not include the individual's current sibling group if he is the only one in it.
ret["AFS"] = clone(AFS);   
    
    
    // now make some copies of what it looked like before the update, because if we just leave the individual
    // in the same full sibling group, it will be easy to just put these values back where they belong.
    NumericVector LMMFS_orig_col = LMMFS( _, fsp);
    NumericVector PMMFS_orig_col = PMMFS( _, fsp);
    NumericVector KidProngs_orig_col = KidProngs( _, fsp);
    
    // divide out his likelihood contribution to the 
    // the LMMFS and recompute PMMFS and KidProngs for his current sibship
    LMMFS( _, fsp) = LMMFS( _, fsp) / LMMI( _, Ind);
ret["MP_Update_Idx1"]  =  update_marriage_posteriors_in_place(LMMFS, PMMFS, UPG, 9, fsp_vec);
    update_marriage_node_kid_prongs_in_place(FSL, PMMFS, KidProngs, 9, 3, TP, fs_vec);
  
 ret["LMMFS_orig"] =  LMMFS_orig_col;
 ret["LMMFS_test"] =  LMMFS( _, fsp);
 
 ret["PMMFS_orig"] =  PMMFS_orig_col;
 ret["PMMFS_test"] =  PMMFS( _, fsp);
 
 ret["KidProngs_orig"] =  KidProngs_orig_col;
 ret["KidProngs_test"] =  KidProngs( _, fsp);

    // now, we can zoom over all the KidProngs for sibling groups that are part of his Acceptable Sibling Groups 
    //and compute the full conditional likelihood of Ind belonging to each. 
    NumericVector FCLs = kid_prongs_times_ind_likelihoods(FSL, IndLiks( _, Ind), KidProngs, AFS);

ret["FCLs"] = FCLs;

    // here is the individual's genotype likelihood given that he is a singleton:
    double solo_lik = geno_post_c(Gfreqs, IndLiks( _, Ind));
    
ret["solo_lik"] = solo_lik;
    
    // here is what we need from our simple pseudo_prior
    List PseudoPri = pseudo_prior(FSL, IFS[Ind], AFS);
    
    // Now we have all the necessary ingredients to compute the full conditional.
    NumericVector PriTimesLik = as<NumericVector>(PseudoPri["afs"]) * FCLs;
    double normo = sum(PriTimesLik);
    double solo_pri = as<double>(PseudoPri["solo"]);
    
    double solo_prob = (solo_pri * solo_lik) /  ( (solo_pri * solo_lik) + ((1.0 - solo_pri) * normo) );
    
    
    if(R::runif(0, 1) < solo_prob) {
 //     ret[6] = "Solo";   // Here we need to fill in what we do when we make their own sibship
      if(n==1) {  // if he was in a singleton sibship to begin with, just put him back there.
        // since we never moved him out of the FSL, we just put the probs back correctly.
        LMMFS( _, fsp) = LMMFS_orig_col;
        PMMFS( _, fsp) = PMMFS_orig_col;
        KidProngs( _, fsp) = KidProngs_orig_col;
        New_fs = fs;
        New_fs_vec = IntegerVector::create(New_fs);
        New_fsp = fsp;
        New_fsp_vec = IntegerVector::create(New_fsp);
      }
      else {
        New_fs = Pile.back();  // If he was not a singleton before, we have to grab a new slot for him from the Pile
        New_fs_vec = IntegerVector::create(New_fs);
        Pile.pop_back();
        New_fsp = MatPile.back();
        New_fsp_vec = IntegerVector::create(New_fsp);
        MatPile.pop_back();
        
        // take him out of the old one
        tmp = FSL[fs];
        IntegerVector xx = as<IntegerVector>(tmp["Indivs"]);
        int xxl = xx.length();
        IntegerVector tmp3(xxl - 1); // length of new vector
        int j=0;
        for(int i=0; i<xxl; i++) {
          if(xx[i] != Ind) tmp3[j++] = xx[i];
        }
        tmp["Indivs"] = tmp3;
        
        // then put him into the FSL that we just pulled off the Pile
        tmp = FSL[New_fs];
        tmp["Indivs"] = NumericVector::create(Ind);
        tmp["LMMI_Idx"] = New_fsp;
        // Then recalculate the probs and things
        LMMFS( _, New_fsp) = LMMI( _, Ind);
        update_marriage_posteriors_in_place(LMMFS, PMMFS, UPG, 9, New_fsp_vec);
        update_marriage_node_kid_prongs_in_place(FSL, PMMFS, KidProngs, 9, 3, TP, New_fs_vec);
        
        // And don't forget to update his IFS entry:
        IFS[Ind] = New_fs;
      }
    }
    else {  // Put Ind into an existing sibship and put his old one on the Pile if n==1
  //    ret[6] = "Siblio";   
      NumericVector probs = PriTimesLik / normo;
      double val = R::runif(0,1);  // the random value
      double sum = 0.0;
      int pn = probs.length();
      int idx = -1;
      for(int i=0; i<pn; i++) {
        sum += probs[i];
        if(sum >= val) {
          idx = i;
          break;
        }
      }
      
      // so, now, AFS[idx] holds the index of the full sibgroup we will assign this guy to, and we might
      New_fs = AFS[idx];


      // now we can deal with recomputing the likelihoods, etc., as need be:
      if(New_fs == fs) { // in this case there was no change and we just copy the values back
        LMMFS( _, fsp) = LMMFS_orig_col;
        PMMFS( _, fsp) = PMMFS_orig_col;
        KidProngs( _, fsp) = KidProngs_orig_col;
      }
      else {  // otherwise there was a change to a new sibgroup.  So, we need to compute the new LMMFS, etc.
            // Note that we already removed him from the likelihood of his previous sibship.
        // First remove him from the old sibship in FSL
        tmp = FSL[fs];
        IntegerVector xx = as<IntegerVector>(tmp["Indivs"]);
        int xxl = xx.length();
        IntegerVector tmp3(xxl - 1); // length of new vector
        int j=0;
        for(int i=0; i<xxl; i++) {
          if(xx[i] != Ind) tmp3[j++] = xx[i];
        }
        tmp["Indivs"] = tmp3;
        
        // Then add him to the vector of Indivs in the new sibgroup
        tmp = FSL[New_fs];
        xx = as<IntegerVector>(tmp["Indivs"]);
        xxl = xx.length();
        IntegerVector tmp2(xxl + 1);  // length of the new vector
        for(int i=0; i<xxl; i++) {
          tmp2[i] = xx[i];
        }
        tmp2[xxl] = Ind;
    //    ret[5] = tmp2;
        tmp["Indivs"] = tmp2;
        
        // then update the probs
        New_fsp = as<int>(tmp["LMMI_Idx"]);
        LMMFS( _, New_fsp) = LMMFS( _, New_fsp) * LMMI( _, Ind);
        update_marriage_posteriors_in_place(LMMFS, PMMFS, UPG, 9, New_fsp_vec);
        update_marriage_node_kid_prongs_in_place(FSL, PMMFS, KidProngs, 9, 3, TP, New_fs_vec);
        
        // And don't forget to update his IFS entry:
        IFS[Ind] = New_fs;
        
        // finally, down here, if n==1, toss his old sibship onto the pile
        if(n==1) {
          Pile.push_back(fs);
          MatPile.push_back(fsp);
        }
      }
    }
    
    
/*    ret[0] = IFS;
    ret[1] = fs;
    ret[2] = fsp;
    ret[3] = New_fs;
    ret[4] = New_fsp;
    ret[7] = Pile;
    ret[8] = MatPile;
    ret[9] = PMMFS;
*/
    ret["Pile"] = Pile;
    ret["MatPile"] = MatPile;
    ret["IFS_end"] = IFS;
    ret["FSL_end"] = FSL;
     
    return(ret);
}




//' return the high-sib-pair-logl indivs for each indiv
//' 
//' I wrote this to see how much faster I could implement the functionality of find_high_logl_sib_pairs
//' using Rcpp.  find_high_logl_sib_pairs is painfully slow, and I think this should be much faster
//' @param FSP a 9 x L matrix (3 x 3 x L) of expected genotype frequencies of a full sibling pair,
//' but we pass it in is a vector.
//' @param UPF a 9 x L matrix of expected genotype frequencies of an unrelated pair, which we also
//' pass in as a vector.
//' @param G an L x N matrix of 0, 1, or 2 or NA giving the genotypes of the individuals 
//' @param loglV the cutoff point above which you will accept the pairs.
//' @export
// [[Rcpp::export]]
List high_logl_pairs(NumericVector FSP, NumericVector UPF, IntegerMatrix G, double loglV) {

  int N = G.ncol();  // number of individuals
  int L = G.nrow();  // number of loci
  NumericVector logFSP = log(FSP);  // we can just log these and add
  NumericVector logUPF = log(UPF);
  List ret(N);  // for returning the values.  
  NumericVector TMP(N); // for storing values or the loglRatio
  double sum;
  int n;
  
  // this is for picking out the elements of FSP and UPF
  #define UF_(g1, g2, l)  g1 + (g2*3) + (9*l)

  // cycle over everyone and do it
  for(int a=0; a<N; a++) {
    n = 0;
    for(int b=0; b<N; b++) {
      sum = 0.0;
      for(int l=0; l<L; l++) {
        if( !(IntegerVector::is_na(G(l, a))) && !(IntegerVector::is_na(G(l, b))) ) {  // if they aren't missing data
          sum += logFSP[UF_(G(l, a), G(l, b), l)] - logUPF[UF_(G(l, a), G(l, b), l)];
        }
      }
      TMP[b] = sum;
      if(TMP[b] > loglV && b != a) n++;
    }
    // now, at this juncture, TMP has all the values for individual a compared to everyone and since n 
    // tells us how many are > loglV, then make an output vector of  length n and
    // copy stuff across.
    IntegerVector reta(n);
    n=0;
    for(int b=0; b<N; b++) {
      if(TMP[b] > loglV && b != a) reta[n++] = b; 
    }
    ret[a] = reta;
  }
  
  return(ret);
  #undef UF_
}











