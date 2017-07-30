#include <algorithm>
#include "Rcpp.h"
#include "UTestFactory.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"
#include "DayBlock.h"

using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

double UTestFactory::epsilon = 3.0;


UTestFactory::UTestFactory(Rcpp::NumericMatrix U,
			   Rcpp::IntegerVector X_rcpp,
			   Rcpp::List w_day_blocks,
			   Rcpp::IntegerVector w_to_days_idx,
			   Rcpp::IntegerVector w_cyc_to_subj_idx,
			   Rcpp::List subj_day_blocks,
			   Rcpp::IntegerVector day_to_subj_idx,
			   Rcpp::List gamma_specs,
			   Rcpp::NumericVector phi_specs,
			   int fw_len,
			   int n_burn,
			   int n_samp,
			   Rcpp::List test_data) :
    // usual input
    U(U),
    X_rcpp(X_rcpp),
    preg_cyc(w_day_blocks),
    w_to_days_idx(w_to_days_idx),
    w_cyc_to_subj_idx(w_cyc_to_subj_idx),
    subj_day_blocks(subj_day_blocks),
    day_to_subj_idx(day_to_subj_idx),
    gamma_specs(gamma_specs),
    phi_specs(phi_specs),
    fw_len(fw_len),
    n_burn(n_burn),
    n_samp(n_samp),
    // testing data
    input_ubeta(as<NumericVector>(test_data["input_ubeta"])),
    input_w(as<IntegerVector>(test_data["input_w"])),
    input_xi(as<NumericVector>(test_data["input_xi"])),
    target_data_gam_cat(as<NumericVector>(test_data["target_data_gam_cat"])),
    target_data_phi(as<NumericVector>(test_data["target_data_phi"])),
    target_samples_gam_cat(as<Rcpp::List>(test_data["target_samples_gam_cat"])),
    target_samples_phi(as<NumericVector>(test_data["target_samples_phi"])),
    target_samples_w(as<NumericVector>(test_data["target_samples_w"])),
    target_samples_xi(as<NumericVector>(test_data["target_samples_xi"])),
    // global testing objects
    seed_vals(as<IntegerVector>(test_data["seed_vals"])),
    // derived data
    n_days(X_rcpp.size()),
    n_subj(subj_day_blocks.size()) {
}


XiGen* UTestFactory::xi() {
    XiGen* xi = new XiGen(subj_day_blocks, n_samp, true);
    std::copy(input_xi.begin(), input_xi.end(), xi->m_vals);
    return xi;
}


XiGen* UTestFactory::xi_no_rec() {
    XiGen* xi_no_rec = new XiGen(subj_day_blocks, n_samp, false);
    std::copy(input_xi.begin(), input_xi.end(), xi_no_rec->m_vals);
    return xi_no_rec;
}


WGen* UTestFactory::W() {
    WGen* W = new WGen(preg_cyc, w_to_days_idx, w_cyc_to_subj_idx, fw_len);
    std::copy(input_w.begin(), input_w.end(), W->m_vals);
    // calculate sums for pregnancy cycles
    int w_ctr = 0;
    int* w_vals = W->m_vals;
    int* w_sums = W->m_sums;
    for (const PregCyc* curr = W->m_preg_cyc; curr < W->m_preg_cyc + W->m_n_preg_cyc; ++curr) {
    	*w_sums = std::accumulate(w_vals + w_ctr, w_vals + w_ctr + curr->n_days, 0.0);
    	++w_sums;
    	w_ctr += curr->n_days;
    }
    return W;
}


PhiGen* UTestFactory::phi() {
    return new PhiGen(phi_specs, n_samp, true);
}


PhiGen* UTestFactory::phi_no_rec() {
    return new PhiGen(phi_specs, n_samp, false);
}


UProdBeta* UTestFactory::ubeta() {
    UProdBeta* ubeta = new UProdBeta(n_days);
    std::copy(input_ubeta.begin(), input_ubeta.end(), ubeta->m_vals);
    ubeta->update_exp(X_rcpp.begin());
    return ubeta;
}


bool UTestFactory::eq_dbl(double a, double b) {
    return (-epsilon < a - b) && (a - b < epsilon);
}
