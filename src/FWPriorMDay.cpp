#include <vector>
#include "Rmath.h"

#include "CoefGen.h"
#include "GammaFWDay.h"
#include "FWPriors.h"

class GammaFWDay;

int sample_multi_index(std::vector<double> probs);




MDay::MDay(Rcpp::NumericVector log_prior_probs_rcpp, int n_samp, bool record_status):
    m_vals_rcpp            {(Rcpp::NumericVector(Rcpp::no_init(record_status ? n_samp : 1)))},
    m_vals                 {(m_vals_rcpp.begin())},
    m_log_prior_probs_rcpp {log_prior_probs_rcpp},
    m_log_prior_probs      {(m_log_prior_probs_rcpp.begin())},
    m_n_days_fw            {static_cast<int>(log_prior_probs_rcpp.size())},
    m_record_status        {record_status}
{}




void MDay::sample(const CoefGen& coefs,
                  const Mu& mu,
                  const Nu& nu) {

    std::vector<double> mday_log_liks(m_n_days_fw);
    std::vector<double> conditional_probs(m_n_days_fw);

    // the maximum observed log-likelihood value
    double max_mday_log_lik = std::numeric_limits<double>::lowest();

    // calculate the log likelihood for each of the possible max-day positions
    // and store the maximum value in `max_mday_log_lik`
    for (int i = 0; i < m_n_days_fw; ++i) {

        mday_log_liks[i] = calc_mday_log_lik(i, coefs, mu, nu) + m_log_prior_probs[i];

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
                               const Nu& nu) {

    // extract priors for convenience
    const double mu_val   = mu.val();
    const double nu_val   = nu.val();

    // each iteration adds in the log-likelihood value for the current FW day
    // coefficient to the running total stored in `sum_log_lik`
    double sum_log_lik = 0.0;
    int k = 0;
    for (GammaGen** curr_gamma = coefs.m_fw_coef_start; curr_gamma < coefs.m_fw_coef_end; ++curr_gamma) {

        // current FW day index and coefficient value
        double curr_gam_val = (*curr_gamma)->m_gam_val;
        double decay_val    = dynamic_cast<GammaFWDay*>(coefs.m_gamma[k++])->decay(proposal_peak_idx);

        sum_log_lik += R::dgamma(curr_gam_val, nu_val, decay_val * mu_val / nu_val, 1);
    }

    return sum_log_lik;
}
