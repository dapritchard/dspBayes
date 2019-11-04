#include "Rcpp.h"
#include "CoefGen.h"
#include "GammaGen.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "FWPriors.h"


#define GAMMA_GEN_TYPE_FW_DAY 3
extern bool g_record_status;




CoefGen::CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp) :
    m_gamma         {GammaGen::create_arr(U, gamma_specs)},
    m_fw_coef_start {nullptr},
    m_fw_coef_end   {nullptr},
    m_vals_rcpp     {Rcpp::NumericVector(Rcpp::no_init(gamma_specs.size() * n_samp))},
    m_vals          {m_vals_rcpp.begin()},
    m_n_psi         {0},
    m_n_gamma       {static_cast<int>(gamma_specs.size())}
{
    int t;

    // find the location of the first fertile window day.  If there are no
    // fertile window days then it takes the value one past the last gamma
    // coefficient.
    for (t = 0; t < gamma_specs.size(); ++t) {
        const Rcpp::NumericVector& curr_gamma_specs = Rcpp::as<Rcpp::NumericVector>(gamma_specs[t]);
        if (curr_gamma_specs["type"] == GAMMA_GEN_TYPE_FW_DAY) {
            break;
        }
    }
    m_fw_coef_start = m_gamma + t;

    // find the location one past the last fertile window day.  If there are no
    // fertile window days then it takes the value one past the last gamma
    // coefficient.
    if (t < gamma_specs.size()) ++t;
    for (; t < gamma_specs.size(); ++t) {
        const Rcpp::NumericVector& curr_gamma_specs = Rcpp::as<Rcpp::NumericVector>(gamma_specs[t]);
        if (curr_gamma_specs["type"] != GAMMA_GEN_TYPE_FW_DAY) {
            break;
        }
    }
    m_fw_coef_end = m_gamma + t;
}




CoefGen::~CoefGen() {
    for (int h = 0; h < m_n_gamma; ++h) {
        delete m_gamma[h];
    }
    delete[] m_gamma;
}




void CoefGen::sample(const WGen& W,
                     const XiGen& xi,
                     UProdBeta& ubeta,
                     const int* X,
                     const FWPriors& fw_priors) {

    // if we're past the burn-in phase then update `m_vals` so that we don't
    // overwrite the previous samples in the current scan
    if (g_record_status) {
        m_vals += m_n_gamma;
    }

    // each iteration updates one gamma_h term and correspondingly adjusts
    // the value of `ubeta`.
    for (int j = 0; j < m_n_gamma; ++j) {
        m_vals[j] = m_gamma[j]->sample(W, xi, ubeta, X, fw_priors);

        // switch (j) {
        // case  0: m_gamma[ 0]->m_gam_val = 0.11; m_gamma[ 0]->m_beta_val = std::log(0.11); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.11)); break;
        // case  1: m_gamma[ 1]->m_gam_val = 0.22; m_gamma[ 1]->m_beta_val = std::log(0.22); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.22)); break;
        // case  2: m_gamma[ 2]->m_gam_val = 0.44; m_gamma[ 2]->m_beta_val = std::log(0.44); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.44)); break;
        // case  3: m_gamma[ 3]->m_gam_val = 0.22; m_gamma[ 3]->m_beta_val = std::log(0.22); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.22)); break;
        // case  4: m_gamma[ 4]->m_gam_val = 0.11; m_gamma[ 4]->m_beta_val = std::log(0.11); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.11)); break;
        // case  5: m_gamma[ 5]->m_gam_val = 0.92; m_gamma[ 5]->m_beta_val = std::log(0.92); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.92)); break;
        // case  6: m_gamma[ 6]->m_gam_val = 0.65; m_gamma[ 6]->m_beta_val = std::log(0.65); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.65)); break;
        // case  7: m_gamma[ 7]->m_gam_val = 0.36; m_gamma[ 7]->m_beta_val = std::log(0.36); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.36)); break;
        // case  8: m_gamma[ 8]->m_gam_val = 0.80; m_gamma[ 8]->m_beta_val = std::log(0.80); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.80)); break;
        // case  9: m_gamma[ 9]->m_gam_val = 0.63; m_gamma[ 9]->m_beta_val = std::log(0.63); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(0.63)); break;
        // case 10: m_gamma[10]->m_gam_val = 2.23; m_gamma[10]->m_beta_val = std::log(2.23); ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh, std::log(2.23)); break;
        // case 11: m_gamma[11]->m_gam_val = 1.00; m_gamma[11]->m_beta_val =            0.0; ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh,            0.0); break;
        // case 12: m_gamma[12]->m_gam_val = 1.00; m_gamma[12]->m_beta_val =            0.0; ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh,            0.0); break;
        // case 13: m_gamma[13]->m_gam_val = 1.00; m_gamma[13]->m_beta_val =            0.0; ubeta.add_uh_prod_beta_h(m_gamma[j]->m_Uh,            0.0); break;
        // }
    }

}
