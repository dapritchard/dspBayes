#include "Rcpp.h"
#include "GammaGen.h"

#define GAMMA_GEN_TYPE_CATEG  0
#define GAMMA_GEN_TYPE_ADAPT  1
#define GAMMA_GEN_TYPE_TUNE   2
#define GAMMA_GEN_TYPE_SMOOTH 3


GammaGen::GammaGen(Rcpp::NumericMatrix& U, Rcpp::NumericVector& coef_specs) :
    // initialization list
    m_beta_val(0),
    m_gam_val(1),
    m_hyp_p(coef_specs["hyp_p"]),
    m_hyp_a(coef_specs["hyp_a"]),
    m_hyp_b(coef_specs["hyp_b"]),
    m_bnd_l(coef_specs["bnd_l"]),
    m_bnd_u(coef_specs["bnd_u"])
    m_Uh(U.begin() + coef_specs("h") * U.nrow()),
    m_n_days(U.nrow())  {
}




GammaGen** GammaGen::list_to_arr(Rcpp::List& gamma_specs) {

    GammaGen** gamma = new GammaGen*[gamma_specs.size()];

    for (int t = 0; t < gamma_specs.size(); ++t) {

	// initialize the appropriate subclass of `GammaGen`, as specified by
	// `gamma_specs[t]`
	switch(gamma_specs[t]["type"]) {
	case GAMMA_GEN_TYPE_CATEG:
	    gamma[t] = new GammaCateg(gamma_specs[t]);
	    break;
	case GAMMA_GEN_TYPE_ADAPT:
	    // ******************** TODO
	    break;
	case GAMMA_GEN_TYPE_TUNE:
	    // ******************** TODO
	    break;
	case GAMMA_GEN_TYPE_SMOOTH:
	    // ******************** TODO
	    break;
	}
    }

    return gamma;
}
