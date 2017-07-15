#include "Rcpp.h"
#include "GammaGen.h"

GammaGen::GammaGen(Rcpp::NumericMatrix& U,
		   Rcpp::init) :
    // initialization list
    m_beta_val(0),
    m_gam_val(1),
    m_hyp_p(init["hyp_p"]),
    m_hyp_a(init["hyp_a"]),
    m_hyp_b(init["hyp_b"]),
    m_bnd_l(init["bnd_l"]),
    m_bnd_u(init["bnd_u"])
    m_Uh(U.begin() + U.nrow() * init("h")),
    m_n_days(U.nrow())  {
}
