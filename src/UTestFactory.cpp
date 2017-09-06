#include <algorithm>
#include "Rcpp.h"

#include "CoefGen.h"
#include "UTestFactory.h"
#include "UGenVar.h"
#include "XGen.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"
#include "DayBlock.h"

using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

double UTestFactory::epsilon = 0.0;


UTestFactory::UTestFactory(Rcpp::NumericMatrix u_rcpp,
			   Rcpp::IntegerVector x_rcpp,
			   Rcpp::List          w_day_blocks,
			   Rcpp::IntegerVector w_to_days_idx,
			   Rcpp::IntegerVector w_cyc_to_subj_idx,
			   Rcpp::List          subj_day_blocks,
			   Rcpp::IntegerVector day_to_subj_idx,
			   Rcpp::List          gamma_specs,
			   Rcpp::NumericVector phi_specs,
			   Rcpp::List          x_miss_cyc,
			   Rcpp::List          x_miss_day,
			   Rcpp::NumericVector utau_rcpp,
			   Rcpp::List          tau_coefs,
			   Rcpp::List          u_miss_info,
			   Rcpp::IntegerVector u_miss_type,
			   Rcpp::IntegerVector u_preg_map,
			   Rcpp::IntegerVector u_sex_map,
			   int fw_len,
			   int n_burn,
			   int n_samp,
			   Rcpp::List test_data) :
    // usual input
    u_rcpp(u_rcpp),
    x_rcpp(x_rcpp),
    preg_cyc(w_day_blocks),
    w_to_days_idx(w_to_days_idx),
    w_cyc_to_subj_idx(w_cyc_to_subj_idx),
    subj_day_blocks(subj_day_blocks),
    day_to_subj_idx(day_to_subj_idx),
    gamma_specs(gamma_specs),
    phi_specs(phi_specs),
    x_miss_cyc(x_miss_cyc),
    x_miss_day(x_miss_day),
    utau_rcpp(utau_rcpp),
    tau_coefs(tau_coefs),
    u_miss_info(u_miss_info),
    u_miss_type(u_miss_type),
    u_preg_map(u_preg_map),
    u_sex_map(u_sex_map),
    fw_len(fw_len),
    n_burn(n_burn),
    n_samp(n_samp),
    // testing data
    input_beta_coefs(as<NumericVector>(test_data["input_beta_coefs"])),
    input_gamma_specs(as<Rcpp::List>(test_data["input_gamma_specs"])),
    input_u_categ(as<Rcpp::List>(test_data["input_u_categ"])),
    input_ubeta(as<NumericVector>(test_data["input_ubeta"])),
    input_w(as<IntegerVector>(test_data["input_w"])),
    input_x(as<NumericVector>(test_data["input_x"])),
    input_xi(as<NumericVector>(test_data["input_xi"])),
    target_data_gamma_categ(as<NumericVector>(test_data["target_data_gamma_categ"])),
    target_data_phi(as<NumericVector>(test_data["target_data_phi"])),
    target_data_u_categ(as<IntegerVector>(test_data["target_data_u_categ"])),
    target_samples_gamma_categ(as<Rcpp::List>(test_data["target_samples_gamma_categ"])),
    target_samples_u_categ(as<Rcpp::List>(test_data["target_samples_u_categ"])),
    target_samples_phi(as<NumericVector>(test_data["target_samples_phi"])),
    target_samples_w(as<NumericVector>(test_data["target_samples_w"])),
    target_samples_x(as<Rcpp::List>(test_data["target_samples_x"])),
    target_samples_xi(as<NumericVector>(test_data["target_samples_xi"])),
    // global testing objects
    seed_vals(as<IntegerVector>(test_data["seed_vals"])),
    // derived data
    n_days(x_rcpp.size()),
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
    ubeta->update_exp(x_rcpp.begin());
    return ubeta;
}


UProdTau* UTestFactory::utau() {
    return new UProdTau(utau_rcpp, tau_coefs);
}


GammaCateg* UTestFactory::gamma_categ_all() {
    GammaCateg* all = new GammaCateg(u_rcpp, as<NumericVector>(input_gamma_specs["gamma_specs_all"]));
    all->m_beta_val = target_data_gamma_categ["beta_prev"];
    return all;
}


GammaCateg* UTestFactory::gamma_categ_zero_one() {
    GammaCateg* zero_one = new GammaCateg(u_rcpp, as<NumericVector>(input_gamma_specs["gamma_specs_zero_one"]));
    zero_one->m_beta_val = target_data_gamma_categ["beta_prev"];
    return zero_one;
}


GammaCateg* UTestFactory::gamma_categ_one_inf() {
    GammaCateg* one_inf = new GammaCateg(u_rcpp, as<NumericVector>(input_gamma_specs["gamma_specs_one_inf"]));
    one_inf->m_beta_val = target_data_gamma_categ["beta_prev"];
    return one_inf;
}


GammaCateg* UTestFactory::gamma_categ_zero_half() {
    GammaCateg* zero_half = new GammaCateg(u_rcpp, as<NumericVector>(input_gamma_specs["gamma_specs_zero_half"]));
    zero_half->m_beta_val = target_data_gamma_categ["beta_prev"];
    return zero_half;
}


XGen* UTestFactory::X() {
    return new XGen(x_rcpp, x_miss_cyc, x_miss_day, tau_coefs["cohort_sex_prob"], tau_coefs["sex_coef"]);
}


int** UTestFactory::X_temp() {
    int** X = new int*;
    *X = x_rcpp.begin();
    return X;
}


UGenVarCateg* UTestFactory::u_categ(Rcpp::NumericMatrix* u_rcpp_copy) {

    Rcpp::List curr_var(u_miss_info[target_data_u_categ["var_idx"]]);

    Rcpp::IntegerVector var_info     ( as<Rcpp::IntegerVector>(curr_var["var_info"])      );
    Rcpp::NumericVector u_prior_probs( as<Rcpp::NumericVector>(curr_var["u_prior_probs"]) );
    Rcpp::List var_block_list        ( as<Rcpp::List>(curr_var["var_block_list"])         );

    return new UGenVarCateg(*u_rcpp_copy, var_info, u_prior_probs, var_block_list, u_preg_map, u_sex_map);
}


CoefGen* UTestFactory::coefs() {

    CoefGen* coefs = new CoefGen(u_rcpp, gamma_specs, n_samp);
    std::copy(input_beta_coefs.begin(), input_beta_coefs.end(), coefs->m_vals);

    return coefs;
}


XGen::XMissDay** UTestFactory::XMissDay() {
    XGen::XMissDay** miss_day = new XGen::XMissDay*;
    *miss_day = XGen::XMissDay::list_to_arr(x_miss_day);
    return miss_day;
}


bool UTestFactory::eq_dbl(double a, double b) {
    return (-epsilon < a - b) && (a - b < epsilon);
}
