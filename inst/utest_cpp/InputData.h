class InputData {

public:

    Rcpp::NumericMatrix U;
    Rcpp::IntegerVector X_rcpp;
    Rcpp::List preg_cyc;
    Rcpp::IntegerVector w_days_idx;
    Rcpp::IntegerVector w_cyc_idx;
    Rcpp::List subj_days;
    Rcpp::IntegerVector subj_idx;
    Rcpp::List gamma_specs;
    Rcpp::NumericVector phi_hyper;
    int fw_len;
    int n_burn;
    int n_samp;

};
