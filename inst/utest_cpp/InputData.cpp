#include <iostream>
#include <string>
#include <algorithm>
#include "Rcpp.h"
#include "InputData.h"
#include "read_data.h"

// the number of values in each list element in the data on file.  For example,
// if there where 4 values in each of 5 list elements, then the data would have
// 20 values.
#define DSP_BAYES_N_DAY_BLOCK_ELEM 4
#define DSP_BAYES_N_SUBJ_BLOCK_ELEM 2
#define DSP_BAYES_N_GAMMA_ELEM 7

using std::cout;
using std::max;
using Rcpp::as;




InputData::InputData(std::string dir_nm) {

    int n_elem;

    int* ptr_int;
    double* ptr_dbl;

    // read `X_rcpp`
    Rcpp::NumericVector X_rcpp = file_to_rcpp_numer(dir_nm.append("/input_X_rcpp"));

    // read `U`
    Rcpp::NumericMatrix U = file_to_rcpp_numer_mat(dir_nm.append("/input_U"), X_rcpp.size());

    // read `w_day_blocks`
    ptr_int = read_data<int>(dir_nm.append("/input_w_day_blocks"), &n_elem);
    Rcpp::List w_day_blocks(n_elem);
    for (int k = 0; k < n_elem; ++k) {

    	int s = k * DSP_BAYES_N_DAY_BLOCK_ELEM;

    	w_day_blocks[k] = Rcpp::IntegerVector::create(Rcpp::Named("beg_idx")  = ptr_int[s + 0],
    						      Rcpp::Named("n_days")   = ptr_int[s + 1],
    						      Rcpp::Named("subj_idx") = ptr_int[s + 2]);
    }
    delete[] ptr_int;

    // read `w_to_days_idx`
    Rcpp::IntegerVector w_to_days_idx = file_to_rcpp_int(dir_nm.append("/input_w_to_days_idx"));

    // // read `w_cyc_to_cyc_idx`
    Rcpp::IntegerVector w_cyc_to_cyc_idx = file_to_rcpp_int(dir_nm.append("/input_w_cyc_to_cyc_idx"));

    // read `subj_day_blocks`
    ptr_int = read_data<int>(dir_nm.append("/input_subj_day_blocks"), &n_elem);
    Rcpp::List subj_day_blocks(n_elem);
    for (int k = 0; k < n_elem; ++k) {

    	int s = k * DSP_BAYES_N_SUBJ_BLOCK_ELEM;

    	subj_day_blocks[k] = Rcpp::IntegerVector::create(Rcpp::Named("beg_idx")  = ptr_int[s + 0],
    							 Rcpp::Named("n_days")   = ptr_int[s + 1]);
    }
    delete[] ptr_int;

    // // read `day_to_subj_idx`
    Rcpp::IntegerVector day_to_subj_idx = file_to_rcpp_int(dir_nm.append("/input_day_to_subj_idx"));

    // read `gamma_specs`
    ptr_dbl = read_data<double>(dir_nm.append("/input_gamma_specs"), &n_elem);
    Rcpp::List gamma_specs(n_elem);
    for (int k = 0; k < n_elem; ++k) {

    	int s = k * DSP_BAYES_N_GAMMA_ELEM;

    	gamma_specs[k] = Rcpp::IntegerVector::create(Rcpp::Named("type")  = ptr_dbl[s + 0],
    						     Rcpp::Named("h")     = ptr_dbl[s + 1],
    						     Rcpp::Named("hyp_a") = ptr_dbl[s + 2],
    						     Rcpp::Named("hyp_b") = ptr_dbl[s + 3],
    						     Rcpp::Named("hyp_p") = ptr_dbl[s + 4],
    						     Rcpp::Named("bnd_l") = ptr_dbl[s + 5],
    						     Rcpp::Named("bnd_u") = ptr_dbl[s + 6]);
    }
    delete[] ptr_dbl;

    // read `phi_specs`
    ptr_dbl = read_data<double>(dir_nm.append("/input_phi_specs"), &n_elem);
    Rcpp::List phi_specs(n_elem);
    phi_specs = Rcpp::NumericVector::create(Rcpp::Named("c1")    = ptr_dbl[0],
    					    Rcpp::Named("c2")    = ptr_dbl[1],
    					    Rcpp::Named("delta") = ptr_dbl[2],
    					    Rcpp::Named("mean")  = ptr_dbl[3]);
    delete[] ptr_dbl;

    // read remaining
    ptr_int = read_data<int>(dir_nm.append("/input_misc_specs"), &n_elem);
    fw_len = ptr_int[0];
    n_burn = ptr_int[1];
    n_samp = ptr_int[2];
    delete[] ptr_int;
}




void InputData::print_visual_verif() const {

    // print first 5 x 5 block of `U`
    cout << "U:  nrow = " << U.nrow() << "  " << U.ncol() <<  "\n";
    cout << "--------------------\n";
    for (int i = 0; i < max((int) U.nrow(), 5); ++i) {
    	for (int j = 0; j < max(U.ncol(), 5); ++j) {
    	    cout << U(i, j) << "  ";
    	}
    	cout << "\n";
    }
    cout << "\n\n";

    // print first 5 elements of `X_rcpp`
    cout << "X:  length = " << X_rcpp.size() << "\n";
    cout << "--------------------\n";
    for (int i = 0; i < max((int) X_rcpp.size(), 5); ++i) {
    	cout << X_rcpp(i) << "  ";
    }
    cout << "\n\n\n";

    // print first two elements of `w_day_blocks`
    cout << "w_day_blocks:  length = " << w_day_blocks.size() << "\n";
    cout << "--------------------\n";
    for (int k = 0; k < max((int) w_day_blocks.size(), 2); ++k) {
    	Rcpp::IntegerVector curr = as<Rcpp::IntegerVector>(w_day_blocks[k]);
    	cout << "beg_idx[" << k << "]: "  << ((int) curr["beg_idx"])  << "  "
    	     << "n_days[" << k << "]: "   << ((int) curr["n_days"])   << "  "
    	     << "subj_idx[" << k << "]: " << ((int) curr["subj_idx"]) << "\n";
    }
    cout << "\n\n";

    // print first 5 elements of `w_to_days_idx`
    cout << "w_to_days_idx:  length = " << w_to_days_idx.size() << "\n";
    cout << "--------------------\n";
    for (int i = 0; i < max((int) w_to_days_idx.size(), 5); ++i) {
    	cout << w_to_days_idx(i) << "  ";
    }
    cout << "\n\n\n";

    // print first 5 elements of `w_cyc_to_cyc_idx`
    cout << "w_cyc_to_cyc_idx:  length = " << w_cyc_to_cyc_idx.size() << "\n";
    cout << "--------------------\n";
    for (int i = 0; i < max((int) w_cyc_to_cyc_idx.size(), 5); ++i) {
    	cout << w_cyc_to_cyc_idx(i) << "  ";
    }
    cout << "\n\n\n";

    // print first two elements of `subj_day_blocks`
    cout << "subj_day_blocks:  length = " << subj_day_blocks.size() << "\n";
    cout << "--------------------\n";
    for (int k = 0; k < max((int) subj_day_blocks.size(), 2); ++k) {
    	Rcpp::IntegerVector curr = as<Rcpp::IntegerVector>(subj_day_blocks[k]);
    	cout << "beg_idx[" << k << "]: " << ((int) curr["beg_idx"]) << "  "
    	     << "n_days[" << k << "]: "  << ((int) curr["n_days"])  << "\n";
    }
    cout << "\n\n";

    // print first 5 elements of `day_to_subj_idx`
    cout << "day_to_subj_idx:  length = " << day_to_subj_idx.size() << "\n";
    cout << "--------------------\n";
    for (int i = 0; i < max((int) day_to_subj_idx.size(), 5); ++i) {
    	cout << day_to_subj_idx(i) << "  ";
    }
    cout << "\n\n\n";

    // print first two elements of `gamma_specs`
    cout << "gamma_specs:  length = " << gamma_specs.size() << "\n";
    cout << "--------------------\n";
    for (int k = 0; k < max((int) gamma_specs.size(), 2); ++k) {
    	Rcpp::IntegerVector curr = as<Rcpp::IntegerVector>(gamma_specs[k]);
    	cout << "type[" << k << "]: "  << ((int) curr["type"])  << "  "
    	     << "h[" << k << "]: "     << ((int) curr["h"])     << "  "
    	     << "hyp_a[" << k << "]: " << ((int) curr["hyp_a"]) << "  "
    	     << "hyp_b[" << k << "]: " << ((int) curr["hyp_b"]) << "  "
    	     << "hyp_p[" << k << "]: " << ((int) curr["hyp_b"]) << "  "
    	     << "bnd_l[" << k << "]: " << ((int) curr["bnd_l"]) << "  "
    	     << "bnd_u[" << k << "]: " << ((int) curr["bnd_u"]) << "\n";
    }
    cout << "\n\n";

    // print `phi_specs`
    cout << "phi_specs:\n";
    cout << "--------------------\n"
    	 << "c1: "    << ((int) phi_specs["c1"])    << "  "
    	 << "c2: "    << ((int) phi_specs["c2"])    << "  "
    	 << "delta: " << ((int) phi_specs["delta"]) << "  "
    	 << "mean: "  << ((int) phi_specs["mean"])  << "\n\n\n";

    // print remaining
    cout << "remaining:\n"
    	 << "--------------------\n"
    	 << "fw_len: " << fw_len << "  "
    	 << "n_burn: " << n_burn << "  "
    	 << "n_samp: " << n_samp << "  ";
}




Rcpp::IntegerVector InputData::file_to_rcpp_int(std::string path) {

    int n;
    int* ptr_start;
    Rcpp::IntegerVector vec;

    ptr_start = read_data<int>(path, &n);
    vec = Rcpp::IntegerVector(ptr_start, ptr_start + n);
    delete[] ptr_start;

    return vec;
}




Rcpp::NumericVector InputData::file_to_rcpp_numer(std::string path) {

    int n;
    double* ptr_start;
    Rcpp::NumericVector vec;

    ptr_start = read_data<double>(path, &n);
    vec = Rcpp::NumericVector(ptr_start, ptr_start + n);
    delete[] ptr_start;

    return vec;
}




Rcpp::NumericMatrix InputData::file_to_rcpp_numer_mat(std::string path, int n_row) {

    int n;
    double* ptr_start;
    Rcpp::NumericMatrix mat;

    ptr_start = read_data<double>(path, &n);
    mat = Rcpp::NumericMatrix(n_row, n / n_row, ptr_start);
    delete[] ptr_start;

    return mat;
}
