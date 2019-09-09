#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "XGen.h"
#include "DayBlock.h"
#include "UProdBeta.h"




XiGen::XiGen(Rcpp::List& subj_day_blocks, int n_samp, bool record_status) :
    // m_vals_rcpp(Rcpp::NumericVector(Rcpp::no_init(subj_day_blocks.size() * (record_status ? n_samp : 1)))),
    m_vals_rcpp(Rcpp::NumericVector(subj_day_blocks.size() * (record_status ? n_samp : 1))),
    m_vals(m_vals_rcpp.begin()),
    m_subj(DayBlock::list_to_arr(subj_day_blocks)),
    m_n_subj(subj_day_blocks.size()),
    m_record_status(record_status)
{
    // initialize values for all subjects to 1 (i.e. no fecundability effect)
    for (int i = 0; i < m_n_subj; ++i) {
        m_vals[i] = 1;
    }
}




XiGen::~XiGen() {
    delete[] m_subj;
}




void XiGen::sample(const WGen& W, const PhiGen& phi, const UProdBeta& ubeta, const XGen& X) {

    const int* w_subj_idx = W.subj_idx();
    const int* w_sum_vals = W.sum_vals();
    const double phi_val = phi.val();
    const int* x_vals = X.vals();
    const double* ubeta_exp_vals = ubeta.exp_vals();

    // if we are past the burn-in phase then move the pointer past the samples
    // so that we don't overwrite them
    if (m_record_status && g_record_status) {
        m_vals += m_n_subj;
    }

    // each iteration samples the i-th value of `xi_i` and stores it `m_xi_vals`
    for (int i = 0; i < m_n_subj; ++i) {

        // critical: remove this prior to final version
        m_vals[i] = 1;
        continue;

        int curr_idx, curr_end;
        double curr_w_sum, curr_sum_exp_ubeta;

        // obtain `sum_jk W_ijk`
        if (i == *w_subj_idx) {
            curr_w_sum = *w_sum_vals++;
            ++w_subj_idx;
        }
        else {
            curr_w_sum = 0;
        }

        // index in the day-specific data of the first day and one past the last
        // day for the current subject
        curr_idx = m_subj[i].beg_idx;
        curr_end = curr_idx + m_subj[i].n_days;

        // obtain `sum_jk { X_ijk * exp( u_{ijk}^T beta ) }`
        curr_sum_exp_ubeta = 0;
        for ( ; curr_idx < curr_end; ++curr_idx) {
            if (x_vals[curr_idx]) {
                curr_sum_exp_ubeta += ubeta_exp_vals[curr_idx];
            }
        }

        // sample new value of `xi_i`
        m_vals[i] = R::rgamma(phi_val + curr_w_sum, 1 / (phi_val + curr_sum_exp_ubeta));
    }
}
