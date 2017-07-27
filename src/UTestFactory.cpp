#include <algorithm>
#include "Rcpp.h"
#include "UTestFactory.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"


UTestFactory::UTestFactory(Rcpp::NumericMatrix U,
			   Rcpp::IntegerVector X_rcpp,
			   Rcpp::List w_day_blocks,
			   Rcpp::IntegerVector w_to_days_idx,
			   Rcpp::IntegerVector w_cyc_to_cyc_idx,
			   Rcpp::List subj_day_blocks,
			   Rcpp::IntegerVector day_to_subj_idx,
			   Rcpp::List gamma_specs,
			   Rcpp::NumericVector phi_specs,
			   int fw_len,
			   int n_burn,
			   int n_samp,
			   Rcpp::NumericVector xi_vals,
			   Rcpp::NumericVector test_data_phi,
			   Rcpp::NumericVector test_data_phi_samples) :
    U(U),
    X_rcpp(X_rcpp),
    preg_cyc(w_day_blocks),
    w_to_days_idx(w_to_days_idx),
    w_cyc_to_cyc_idx(w_cyc_to_cyc_idx),
    subj_day_blocks(subj_day_blocks),
    day_to_subj_idx(day_to_subj_idx),
    gamma_specs(gamma_specs),
    phi_specs(phi_specs),
    fw_len(fw_len),
    n_burn(n_burn),
    n_samp(n_samp),
    xi_vals(xi_vals),
    test_data_phi(test_data_phi),
    test_data_phi_samples(test_data_phi_samples),

    n_days(X_rcpp.size()) {
}



XiGen* UTestFactory::xi() {
    XiGen* xi = new XiGen(subj_day_blocks, n_samp, true);
    // TODO: copy vals
    return xi;
}


XiGen* UTestFactory::xi_no_rec() {
    XiGen* xi_no_rec = new XiGen(subj_day_blocks, n_samp, false);
    // TODO: copy vals
    return xi_no_rec;
}


WGen* UTestFactory::W() {
    WGen* W = new WGen(preg_cyc, w_to_days_idx, w_cyc_to_cyc_idx, fw_len);
    // std::copy(w_vals.begin(), w_vals.end(), W.m_vals);
    return W;
}


PhiGen* UTestFactory::phi() {
    PhiGen* phi = new PhiGen(phi_specs, n_samp, true);
    // TODO: copy vals
    return phi;
}


PhiGen* UTestFactory::phi_no_rec() {
    PhiGen* phi_no_rec = new PhiGen(phi_specs, n_samp, false);
    // TODO: copy vals
    return phi_no_rec;
}


UProdBeta* UTestFactory::ubeta() {
    UProdBeta* ubeta = new UProdBeta(n_days);
    // std::copy(ubeta.begin(), ubeta.end(), ubeta->m_vals);
    return ubeta;
}
