#include <vector>
#include "Rmath.h"

#include "CoefGen.h"
#include "GammaFWDay.h"
#include "FWPriors.h"

class GammaFWDay;

int sample_multi_index(std::vector<double> probs);




MDay::MDay(int n_samp, int n_days_fw, bool record_status) :
    m_vals_rcpp     {(Rcpp::NumericVector(Rcpp::no_init(record_status ? n_samp : 1)))},
    m_vals          {(m_vals_rcpp.begin())},
    m_n_days_fw     {n_days_fw},
    m_record_status {record_status}
{}




void MDay::sample(const CoefGen& coefs,
                  const Mu& mu,
                  const Nu& nu,
                  const Delta& delta) {

    std::vector<double> mday_log_liks(m_n_days_fw);
    std::vector<double> conditional_probs(m_n_days_fw);

    // the maximum observed log-likelihood value
    double max_mday_log_lik = std::numeric_limits<double>::lowest();

    // calculate the log likelihood for each of the possible max-day positions
    // and store the maximum value in `max_mday_log_lik`
    for (int i = 0; i < m_n_days_fw; ++i) {

        mday_log_liks[i] = calc_mday_log_lik(i, coefs, mu, nu, delta);

        if (mday_log_liks[i] > max_mday_log_lik) {
            max_mday_log_lik = mday_log_liks[i];
        }
    }

    // subtract away the maximum value from each log likelihood and add in the
    // the corresponding term to `\sum_{i=1}^K \exp{\log_{x_i} - a}` where `a =
    // \max{\log_{x_1}, ..., \log{x_K}}`
    double sum_exp_mday_log_liks = 0.0;
    for (int i = 0; i < m_n_days_fw; ++i) {
        mday_log_liks[i] -= max_mday_log_lik;
        sum_exp_mday_log_liks += std::exp(mday_log_liks[i]);
    }
    double log_sum_exp_mday_log_liks = std::log(sum_exp_mday_log_liks);

    // calculate the conditional probabilities
    for (int i = 0; i < m_n_days_fw; ++i) {
        conditional_probs[i] = std::exp(mday_log_liks[i] - log_sum_exp_mday_log_liks);
    }

    // if we are past the burn-in phase then move the pointer past the samples
    // so that we don't overwrite them
    if (m_record_status && g_record_status) {
        ++m_vals;
    }

    // sample new value and update
    *m_vals = sample_multi_index(conditional_probs);
}




// sample a multinomial value
int sample_multi_index(std::vector<double> probs) {
    double x = R::unif_rand();
    int i;
    for (i = 0; ; ++i) {
        x -= probs[i];
        if (x < 0) break;
    }
    return i;
}




double MDay::calc_mday_log_lik(int proposal_peak_idx,
                               const CoefGen& coefs,
                               const Mu& mu,
                               const Nu& nu,
                               const Delta& delta) {

    // extract priors for convenience
    const double mu_val    = mu.val();
    const double nu_val    = nu.val();
    const double delta_val = delta.val();

    // each iteration adds in the log-likelihood value for the current FW day
    // coefficient to the running total stored in `sum_log_lik`
    double sum_log_lik = 0.0;
    for (GammaGen** curr_gamma = coefs.m_fw_coef_start; curr_gamma < coefs.m_fw_coef_end; ++curr_gamma) {

        // current FW day index and coefficient value
        int curr_day_idx    = dynamic_cast<GammaFWDay*>(*curr_gamma)->m_day_idx;
        double curr_gam_val = (*curr_gamma)->m_gam_val;

        // calculate `delta^{|k-m|}`
        double day_dist  = abs(curr_day_idx - proposal_peak_idx);
        double delta_pow = pow(delta_val, day_dist);

        sum_log_lik += R::dgamma(curr_gam_val, nu_val, delta_pow * mu_val / nu_val, 1);
    }

    return sum_log_lik;
}
