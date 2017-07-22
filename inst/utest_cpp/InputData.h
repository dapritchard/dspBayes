#include <string>
#include "Rcpp.h"


class InputData {

public:

    Rcpp::NumericMatrix U;
    Rcpp::IntegerVector X_rcpp;
    Rcpp::List w_day_blocks;
    Rcpp::IntegerVector w_to_days_idx;
    Rcpp::IntegerVector w_cyc_to_cyc_idx;
    Rcpp::List subj_day_blocks;
    Rcpp::IntegerVector day_to_subj_idx;
    Rcpp::List gamma_specs;
    Rcpp::NumericVector phi_specs;
    int fw_len;
    int n_burn;
    int n_samp;

    InputData(std::string dir_nm);
    void print_visual_verif() const;

    static Rcpp::IntegerVector file_to_rcpp_int(std::string path);
    static Rcpp::NumericVector file_to_rcpp_numer(std::string path);
    static Rcpp::NumericMatrix file_to_rcpp_numer_mat(std::string path, int n_row);
};
