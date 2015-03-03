// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// per_kid_marriage_likelihoods
NumericVector per_kid_marriage_likelihoods(NumericVector offspring_likelihoods, NumericVector transmission_probs);
RcppExport SEXP fullsniplings_per_kid_marriage_likelihoods(SEXP offspring_likelihoodsSEXP, SEXP transmission_probsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type offspring_likelihoods(offspring_likelihoodsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type transmission_probs(transmission_probsSEXP );
        NumericVector __result = per_kid_marriage_likelihoods(offspring_likelihoods, transmission_probs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// C_full_sibling_pair_gfreqs
NumericVector C_full_sibling_pair_gfreqs(int L, int G, NumericVector PP, NumericVector TP, NumericVector GE);
RcppExport SEXP fullsniplings_C_full_sibling_pair_gfreqs(SEXP LSEXP, SEXP GSEXP, SEXP PPSEXP, SEXP TPSEXP, SEXP GESEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type L(LSEXP );
        Rcpp::traits::input_parameter< int >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type PP(PPSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type TP(TPSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type GE(GESEXP );
        NumericVector __result = C_full_sibling_pair_gfreqs(L, G, PP, TP, GE);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// multi_kid_marriage_likelihoods
NumericMatrix multi_kid_marriage_likelihoods(List S, NumericMatrix PK);
RcppExport SEXP fullsniplings_multi_kid_marriage_likelihoods(SEXP SSEXP, SEXP PKSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type S(SSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type PK(PKSEXP );
        NumericMatrix __result = multi_kid_marriage_likelihoods(S, PK);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// update_marriage_likelihoods_in_place
void update_marriage_likelihoods_in_place(List S, NumericMatrix PK, NumericMatrix ML, IntegerVector bz_idx);
RcppExport SEXP fullsniplings_update_marriage_likelihoods_in_place(SEXP SSEXP, SEXP PKSEXP, SEXP MLSEXP, SEXP bz_idxSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type S(SSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type PK(PKSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type ML(MLSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type bz_idx(bz_idxSEXP );
        update_marriage_likelihoods_in_place(S, PK, ML, bz_idx);
    }
    return R_NilValue;
END_RCPP
}
// update_marriage_posteriors_in_place
int update_marriage_posteriors_in_place(NumericMatrix ML, NumericMatrix MP, NumericMatrix Pri, int NGS, IntegerVector bz_idx);
RcppExport SEXP fullsniplings_update_marriage_posteriors_in_place(SEXP MLSEXP, SEXP MPSEXP, SEXP PriSEXP, SEXP NGSSEXP, SEXP bz_idxSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type ML(MLSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type MP(MPSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type Pri(PriSEXP );
        Rcpp::traits::input_parameter< int >::type NGS(NGSSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type bz_idx(bz_idxSEXP );
        int __result = update_marriage_posteriors_in_place(ML, MP, Pri, NGS, bz_idx);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// update_marriage_node_kid_prongs_in_place
void update_marriage_node_kid_prongs_in_place(List S, NumericMatrix MP, NumericMatrix KP, int NGS_P, int NGS_K, NumericVector TP, IntegerVector bz_idx);
RcppExport SEXP fullsniplings_update_marriage_node_kid_prongs_in_place(SEXP SSEXP, SEXP MPSEXP, SEXP KPSEXP, SEXP NGS_PSEXP, SEXP NGS_KSEXP, SEXP TPSEXP, SEXP bz_idxSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type S(SSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type MP(MPSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type KP(KPSEXP );
        Rcpp::traits::input_parameter< int >::type NGS_P(NGS_PSEXP );
        Rcpp::traits::input_parameter< int >::type NGS_K(NGS_KSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type TP(TPSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type bz_idx(bz_idxSEXP );
        update_marriage_node_kid_prongs_in_place(S, MP, KP, NGS_P, NGS_K, TP, bz_idx);
    }
    return R_NilValue;
END_RCPP
}
// kid_prongs_times_ind_likelihoods
NumericVector kid_prongs_times_ind_likelihoods(List FSL, NumericVector IndGenoLik, NumericMatrix KidProngs, IntegerVector AFS);
RcppExport SEXP fullsniplings_kid_prongs_times_ind_likelihoods(SEXP FSLSEXP, SEXP IndGenoLikSEXP, SEXP KidProngsSEXP, SEXP AFSSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type FSL(FSLSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type IndGenoLik(IndGenoLikSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type KidProngs(KidProngsSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type AFS(AFSSEXP );
        NumericVector __result = kid_prongs_times_ind_likelihoods(FSL, IndGenoLik, KidProngs, AFS);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// possible_sibgroups
IntegerVector possible_sibgroups(IntegerVector IFS, IntegerVector AFS);
RcppExport SEXP fullsniplings_possible_sibgroups(SEXP IFSSEXP, SEXP AFSSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type IFS(IFSSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type AFS(AFSSEXP );
        IntegerVector __result = possible_sibgroups(IFS, AFS);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pseudo_prior
List pseudo_prior(List FSL, int IndG, IntegerVector AFS);
RcppExport SEXP fullsniplings_pseudo_prior(SEXP FSLSEXP, SEXP IndGSEXP, SEXP AFSSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type FSL(FSLSEXP );
        Rcpp::traits::input_parameter< int >::type IndG(IndGSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type AFS(AFSSEXP );
        List __result = pseudo_prior(FSL, IndG, AFS);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// geno_post_c
double geno_post_c(NumericVector Gfreqs, NumericVector Liks);
RcppExport SEXP fullsniplings_geno_post_c(SEXP GfreqsSEXP, SEXP LiksSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Gfreqs(GfreqsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Liks(LiksSEXP );
        double __result = geno_post_c(Gfreqs, Liks);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gibbs_update_one_indiv_in_place
List gibbs_update_one_indiv_in_place(List FSL, IntegerVector IFS, NumericMatrix LMMI, NumericMatrix LMMFS, NumericMatrix PMMFS, NumericMatrix KidProngs, std::vector<int> Pile, std::vector<int> MatPile, List AFSL, NumericMatrix Gfreqs, NumericMatrix UPG, NumericVector TP, NumericMatrix IndLiks, int Ind);
RcppExport SEXP fullsniplings_gibbs_update_one_indiv_in_place(SEXP FSLSEXP, SEXP IFSSEXP, SEXP LMMISEXP, SEXP LMMFSSEXP, SEXP PMMFSSEXP, SEXP KidProngsSEXP, SEXP PileSEXP, SEXP MatPileSEXP, SEXP AFSLSEXP, SEXP GfreqsSEXP, SEXP UPGSEXP, SEXP TPSEXP, SEXP IndLiksSEXP, SEXP IndSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type FSL(FSLSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type IFS(IFSSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type LMMI(LMMISEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type LMMFS(LMMFSSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type PMMFS(PMMFSSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type KidProngs(KidProngsSEXP );
        Rcpp::traits::input_parameter< std::vector<int> >::type Pile(PileSEXP );
        Rcpp::traits::input_parameter< std::vector<int> >::type MatPile(MatPileSEXP );
        Rcpp::traits::input_parameter< List >::type AFSL(AFSLSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type Gfreqs(GfreqsSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type UPG(UPGSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type TP(TPSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type IndLiks(IndLiksSEXP );
        Rcpp::traits::input_parameter< int >::type Ind(IndSEXP );
        List __result = gibbs_update_one_indiv_in_place(FSL, IFS, LMMI, LMMFS, PMMFS, KidProngs, Pile, MatPile, AFSL, Gfreqs, UPG, TP, IndLiks, Ind);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// high_logl_pairs
List high_logl_pairs(NumericVector FSP, NumericVector UPF, IntegerMatrix G, double loglV);
RcppExport SEXP fullsniplings_high_logl_pairs(SEXP FSPSEXP, SEXP UPFSEXP, SEXP GSEXP, SEXP loglVSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type FSP(FSPSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type UPF(UPFSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type G(GSEXP );
        Rcpp::traits::input_parameter< double >::type loglV(loglVSEXP );
        List __result = high_logl_pairs(FSP, UPF, G, loglV);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}