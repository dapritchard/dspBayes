#ifndef DSP_BAYES_GAMMA_FWDAY_H_
#define DSP_BAYES_GAMMA_FWDAY_H_

#include "Rcpp.h"
#include "FWPriors.h"


class GammaFWDay : public GammaContMH {
public:

    // // the current value of gamma_h
    // double m_beta_val;
    // double m_gam_val;

    // // the FW day index
    // const int m_day_idx;
    int m_day_idx;

    // the proposal distribution
    double (*m_proposal_fcn)(double cond, double delta);

    GammaFWDay(const Rcpp::NumericMatrix& U,
	       const Rcpp::NumericVector& gamma_specs,
	       int day_idx);
    double sample(const WGen& W,
		  const XiGen& xi,
		  UProdBeta& ubeta,
		  const int* X,
		  const FWPriors& fw_priors);
    double get_log_r(const WGen& W,
		     const XiGen& xi,
		     const UProdBeta& ubeta,
		     const int* X,
		     double proposal_beta,
		     double proposal_gam,
		     const FWPriors& fw_priors);
    double get_gam_log_lik(double proposal_beta,
			   double proposal_gam,
			   FWPriors fw_priors);
};


#endif
