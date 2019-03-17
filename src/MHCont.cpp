#include <cmath>
#include "Rcpp.h"
#include "MHCont.h"


MHCont::MHCont(int n_samp, bool record_status, double proposal_dispersion) :
    m_vals_rcpp     {(Rcpp::NumericVector(Rcpp::no_init(record_status ? n_samp : 1)))},
    m_vals          {(m_vals_rcpp.begin())},
    m_prp_disp      {proposal_dispersion},
    m_accept_ctr    {(0)},
    m_record_status {(record_status)}
{}




// randomly update the parameter value by either accepting the proposal value
// for phi or by keeping the current value.  The proposal value is accepted with
// probability min(1, r).

double MHCont::update(double log_r, double proposal_val) {

    double out;

    // case: accept proposal value
    if ((log_r >= 0) || (std::log(R::unif_rand()) < log_r)) {
        ++m_accept_ctr;
        out = proposal_val;
    }
    // case: reject proposal value
    else {
        out = *m_vals;
    }

    return out;
}
