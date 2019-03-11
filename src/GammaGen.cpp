#include "Rcpp.h"
#include "GammaGen.h"
#include "GammaFWDay.h"

#define GAMMA_GEN_TYPE_CATEG       0
#define GAMMA_GEN_TYPE_CONT_MH     1
#define GAMMA_GEN_TYPE_CONT_ADAPT  2
#define GAMMA_GEN_TYPE_FW_DAY      3

using Rcpp::NumericVector;
using Rcpp::as;




GammaGen::GammaGen(const Rcpp::NumericMatrix& U,
                   const Rcpp::NumericVector& gamma_specs) :
    // initialization list
    m_beta_val(0),
    m_gam_val(1),
    m_hyp_a(gamma_specs["hyp_a"]),
    m_hyp_b(gamma_specs["hyp_b"]),
    m_hyp_p(gamma_specs["hyp_p"]),
    m_bnd_l(gamma_specs["bnd_l"]),
    m_bnd_u(gamma_specs["bnd_u"]),
    m_Uh(U.begin() + ((int) gamma_specs["h"]) * U.nrow()),
    m_n_days(U.nrow())  {
}




GammaGen** GammaGen::create_arr(const Rcpp::NumericMatrix& U,
                                const Rcpp::List& gamma_specs) {

    GammaGen** gamma = new GammaGen*[gamma_specs.size()];

    for (int t = 0; t < gamma_specs.size(); ++t) {

        // convert `List` element to `NumericVector`
        const NumericVector& curr_gamma_specs = as<NumericVector>(gamma_specs[t]);

        // initialize the appropriate subclass of `GammaGen`, as specified by
        // `curr_gamma_specs`
        switch((int) curr_gamma_specs["type"]) {
        case GAMMA_GEN_TYPE_CATEG:
            gamma[t] = new GammaCateg(U, curr_gamma_specs);
            break;
        case GAMMA_GEN_TYPE_CONT_MH:
            gamma[t] = new GammaContMH(U, curr_gamma_specs);
            break;
        case GAMMA_GEN_TYPE_FW_DAY:
            gamma[t] = new GammaFWDay(U, curr_gamma_specs);
            break;
        }
    }

    return gamma;
}
