#include <string>
#include "Rcpp.h"
#include "InputData.h"
#include "read_data.h"

#define DSP_BAYES_N_DAY_BLOCK_ELEM 4
#define DSP_BAYES_N_SUBJ_DAY_BLOCK_ELEM 2
#define DSP_BAYES_N_GAMMA_ELEM 2




InputData::InputData(std::string dir_nm) {

    // int n_days, n_coefs, n_cycs, n_elem;

    // int* ptr_x_rcpp, *ptr_w_day_blocks, *ptr_w_to_days_idx, *ptr_w_cyc_to_cyc_idx,
    // 	*ptr_subj_day_blocks, *ptr_misc_specs;
    // double* ptr_u, *ptr_gamma_specs, *ptr_phi_specs;

    // read `X_rcpp`
    Rcpp::NumericVector X_rcpp = file_to_rcpp_numer(dir_nm.append("/input_X_rcpp"));

    // read `U`
    ptr_u = read_data<double>(dir_nm.append("/input_U"), &n_coefs);
    n_coefs /= n_days;
    U = Rcpp::NumericMatrix(n_days, n_coefs, ptr_u);

    // read `w_day_blocks`
    ptr_w_day_blocks = read_data<int>(dir_nm.append("/input_w_day_blocks"), &n_elem);
    // w_day_blocks = Rcpp::List(n_elem);
    // for (int* curr = ptr_w_day_blocks;
    // 	 curr < ptr_w_day_blocks + n_elem;
    // 	 curr += DSP_BAYES_N_DAY_BLOCK_ELEM) {

    // 	w_day_blocks[i] = Rcpp::IntegerVector::create(Rcpp::Named["beg_idx"]  = curr[0],
    // 						      Rcpp::Named["n_days"]   = curr[1],
    // 						      Rcpp::Named["subj_idx"] = curr[2]);
    // }

    // read `w_to_days_idx`
    Rcpp::IntegerVector w_to_days_idx = file_to_rcpp_int(dir_nm.append("/input_w_to_days_idx"));

    // // read `w_cyc_to_cyc_idx`
    Rcpp::IntegerVector w_cyc_to_cyc_idx = file_to_rcpp_int(dir_nm.append("/input_w_cyc_to_cyc_idx"));

    // // read `subj_day_blocks`
    // ptr_subj_day_blocks = read_data<int>(dir_nm + "/input_subj_day_blocks", &n_elem);
    // subj_day_blocks = Rcpp::List(n_elem);
    // for (int* curr = ptr_w_day_blocks;
    // 	 curr < ptr_w_day_blocks + n_elem;
    // 	 curr +=  DSP_BAYES_N_SUBJ_DAY_BLOCK_ELEM) {

    // 	subj_day_blocks[i] = Rcpp::IntegerVector::create(Rcpp::Named["beg_idx"]  = curr[0],
    // 							 Rcpp::Named["n_days"]   = curr[1]);
    // }

    // // read `day_to_subj_idx`
    Rcpp::IntegerVector day_to_subj_idx = file_to_rcpp_int(dir_nm.append("/input_day_to_subj_idx"));

    // // read `gamma_specs`
    // ptr_gamma_specs = read_data<double>(dir_nm + "/input_gamma_specs", &n_elem);
    // gamma_specs = Rcpp::List(n_elem);
    // for (int* curr = ptr_w_day_blocks;
    // 	 curr < ptr_w_day_blocks + n_elem;
    // 	 curr += DSP_BAYES_N_GAMMA_ELEM) {

    // 	gamma_specs[i] = Rcpp::IntegerVector::create(Rcpp::Named["beg_idx"]  = curr[0],
    // 						     Rcpp::Named["n_days"]   = curr[1],
    // 						     Rcpp::Named["subj_idx"] = curr[2]);
    // }

    // read `phi_specs`
    ptr_phi_specs = read_data<double>(dir_nm + "/input_phi_specs", &n_elem);
    phi_specs = Rcpp::List(n_elem);
    phi_specs = Rcpp::NumericVector::create(Rcpp::Named["c1"]    = ptr_phi_specs[0],
    					    Rcpp::Named["c2"]    = ptr_phi_specs[1],
    					    Rcpp::Named["delta"] = ptr_phi_specs[2],
    					    Rcpp::Named["mean"]  = ptr_phi_specs[3]);

    // // read remaining
    // ptr_misc_specs = read_data<int>(dir_nm + "/input_misc_specs", &n_elem);
    // fw_len = ptr_misc_specs[0];
    // n_burn = ptr_misc_specs[1];
    // n_samp = ptr_misc_specs[2];
}




Rcpp::IntegerVector file_to_rcpp_int(std::string path) {

    int n;
    int* ptr_start;
    Rcpp::IntegerVector vec;

    ptr_start = read_data<int>(path, &n);
    vec = Rcpp::IntegerVector(ptr_start, ptr_start + n);
    delete[] ptr_start;

    return vec;
}




Rcpp::NumericVector file_to_rcpp_numer(std::string path) {

    int n;
    double* ptr_start;
    Rcpp::NumericVector vec;

    ptr_start = read_data<double>(path, &n);
    vec = Rcpp::NumericVector(ptr_start, ptr_start + n);
    delete[] ptr_start;

    return vec;
}
