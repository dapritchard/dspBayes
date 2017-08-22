#include <cmath>
#include <algorithm>
#include "Rcpp.h"

#include "CoefGen.h"
#include "UGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XiGen.h"

#define REF_CELL         -1
#define IN_NON_PREG_CYC  -1
#define NONMISS_SEX_DAY  -1




UGen::UGen() :
    m_col_start(0),
    m_col_end(0),
    m_ref_col(0),
    m_n_categs(0),
    // m_is_ref_cell_coding(0),
    m_miss_block(0)// ,
    // m_w_idx(0)
{
}



// void UGen::sample() {

//     m_curr_block = m_miss_block;
//     m_curr_day = m_miss_day;

//     for ( ; m_curr_block < m_end_block; ++m_curr_block) {

// 	sample_block();
//     }

// }




// void UGen::sample_block() {

//     TODO: check that there are any missing in X

//     // calc p(W | ubeta) perms

//     // calc p(X | utau) perms

//     // calc p(U) perms

//     // sample missing covariate

//     // update miss_block

//     // update ubeta

//     // update utau

//     // -------------------------------

//     double* posterior_w_probs[m_n_categs];
//     double* posterior_x_probs[m_n_categs];

//     calc_posterior_w_probs(posterior_w_probs, W, xi, coefs, ubeta, miss_block);


// }




void UGen::calc_posterior_w(double* posterior_w_probs,
			    double* alt_exp_ubeta_vals,
			    const WGen& W,
			    const XiGen& xi,
			    const CoefGen& coefs,
			    const UProdBeta& ubeta,
			    const UMissBlock* const miss_block) const {

    // const int block_beg_day_idx = miss_block->beg_day_idx;
    // const int block_beg_w_idx   = miss_block->beg_w_idx;
    const int block_n_days      = miss_block->n_days;
    const int block_u_col       = miss_block->u_col;

    const double* block_exp_ubeta_vals = ubeta.exp_vals() + miss_block->beg_day_idx;
    const int* block_w_idx           = m_w_idx + miss_block->beg_w_idx;

    const double* beta_coefs = coefs.vals();
    const int* w_vals        = W.vals();
    const double xi_i        = xi.vals()[miss_block->subj_idx];

    // each iteration calculates `p(W | U, data)` for the current category of
    // the missing covariate
    for (int j = m_col_start; j < m_col_end; ++j) {

	double exp_beta_diff;
	double log_dpois_sum;

	// `exp_beta_diff` is the value of `exp(beta_j - beta^{*})` where
	// `beta^{*}` is the coefficient of beta that corresponds to the
	// category of the missing covariate chosen for the previous sample.
	// Either `beta_j` or `beta^{*}` or both may correspond the the
	// reference cell category, in which case the corresponding value is
	// implicitely 0.
	//
	// case: the j-th column is the one corresponding to the previous choice
	// of category for the missing covariate
	if (j == block_u_col) {
	    exp_beta_diff = 1;
	}
	// case: the j-th column is not the one corresponding to the previous
	// choice of category for the missing covariate
	else {
	    const double beta_j = (j == m_ref_col) ?
		0.0 :
		beta_coefs[j];
	    const double beta_star = (block_u_col == REF_CELL) ?
		0.0 :
		beta_coefs[block_u_col];

	    exp_beta_diff = exp(beta_j - beta_star);
	}

	// each iteration calculated the value of `log p(W_ijk | U, data)` and
	// adds the value to `log_dpois_sum` for `ijk` one of the values of `W`
	// affected by the missing covariate
	for (int r = 0 ; r < block_n_days; ++r) {

	    *alt_exp_ubeta_vals = block_exp_ubeta_vals[r] * exp_beta_diff;

	    // the mean value for the Poisson distribution of `W_ijk`, and the
	    // index in `W.vals()` corresponding to observation `ijk`
	    const double curr_mean_val = xi_i * *alt_exp_ubeta_vals++;
	    const int curr_w_idx       = block_w_idx[r];

	    // the value of `log p(W_ijk | U, data)`.  If `curr_w_idx` has the
	    // flag value indicating that observation `ijk` corresponds to a
	    // day in a non-pregnancy cycle, then `W_ijk` has a value of 0.
	    log_dpois_sum += (curr_w_idx == IN_NON_PREG_CYC) ?
		-curr_mean_val :
		R::dpois(w_vals[curr_w_idx], curr_mean_val, 1);
	}

	// exponentiate to recover the posterior probability, and store the
	// result
	*posterior_w_probs++ = exp(log_dpois_sum);
    }
}




void UGen::calc_posterior_x(double* posterior_x_probs,
			    double* alt_utau_vals,
			    const XGen& X,
			    const UProdTau& utau,
			    const UMissBlock* const miss_block) const {

    const int block_beg_sex_idx = miss_block->beg_sex_idx;
    const int block_n_sex_days  = miss_block->n_sex_days;
    const int block_u_col       = miss_block->u_col;

    const int* x_vals                = X.vals();
    const XGen::XMissDay* x_miss_day = X.miss_day();
    const double sex_coef            = X.sex_coef();

    const double* tau_coefs = utau.coefs();
    const double* utau_vals = utau.vals();

    // each iteration calculates `P(X | U)` corresponding to a value of 1 for
    // the dummy coding in the j-th variable for the observations of `X` that
    // are affected by the missing covariate
    for (int j = m_col_start; j < m_col_end; ++j) {

	double neg_log_probs_sum = 0.0;

	// the value of `tau_j - tau^{*}` where `tau^{*}` is the coefficient of
	// tau that corresponds to the category of `U` chosen for the previous
	// sample.  Either `tau_j` or `tau^{*}` or both may correspond the the
	// reference cell category, in which case the corresponding value is
	// implicitely 0.
	const double tau_j = (j == m_ref_col) ?
	    0.0 :
	    tau_coefs[j];
	const double tau_diff = (block_u_col == REF_CELL) ?
	    tau_j:
	    tau_j - tau_coefs[block_u_col];

	// each iteration calculated the value of `-log P(X_ijk | U_ijk)` and
	// adds the value to `neg_log_probs_sum` for `ijk` one of the values of
	// `X` affected by the missing covariate
	for (int r = 0; r < block_n_sex_days; ++r) {

	    const int curr_x_day_idx = m_x_idx[block_beg_sex_idx + r];
	    const int curr_day_idx   = x_miss_day[curr_x_day_idx].idx;
	    const int curr_sex_prev  = x_miss_day[curr_x_day_idx].prev;

	    // a value of <= 0 corresponds to no sex on the previous day, and >=
	    // 1 corresponds to yes sex on the previous day
	    *alt_utau_vals = utau_vals[curr_x_day_idx] + tau_diff;
	    const double curr_utau_val = (curr_sex_prev <= 0) ?
		*alt_utau_vals++ :
		*alt_utau_vals++ + sex_coef;

	    // add `-log P(X_ijk | U_ijk)` to the running total, where `P(X_ijk
	    // | U_ijk)` follows a logistic regression model
	    neg_log_probs_sum += ((curr_day_idx != IN_NON_PREG_CYC) && x_vals[curr_day_idx]) ?
		log(1 + exp(- curr_utau_val)) :
		log(1 + exp(curr_utau_val));
	}

	// exponentiate the negative of the result to recover the posterior
	// probability, and store the result
	*posterior_x_probs++ = exp(- neg_log_probs_sum);
    }
}




int UGen::sample_covariate(const double* posterior_w_probs,
			   const double* posterior_x_probs) const {

    double unnormalized_probs[m_n_categs];
    double normalizing_constant;

    // each iteration calculates the unnormalized posterior probability for the
    // j-th category for the missing covariate, and adds the value to the
    // running total for the normalizing constant
    normalizing_constant = 0.0;
    for (int j = 0; j < m_n_categs; ++j) {
	normalizing_constant +=
	    unnormalized_probs[j] =
	    posterior_x_probs[j] * posterior_w_probs[j] * m_u_prior_probs[j];
    }

    // sample a value from a `unif(0, normalizing_constant)` distribution
    const double u = R::unif_rand() * normalizing_constant;

    // draw from a multinomial distribution with 1 draw and `m_n_categs` bins
    int j = 0;
    double bin_rhs = unnormalized_probs[0];
    while (u > bin_rhs) {
	bin_rhs += unnormalized_probs[++j];
    }

    return j;
}




void UGen::update_ubeta(UProdBeta& ubeta,
			int u_categ,
			const double* alt_exp_ubeta_vals,
			const UMissBlock* const miss_block) {

    const int block_beg_day_idx = miss_block->beg_day_idx;
    const int block_n_days      = miss_block->n_days;

    double* block_ubeta_vals             = ubeta.vals() + block_beg_day_idx;
    double* block_exp_ubeta_vals         = ubeta.exp_vals() + block_beg_day_idx;
    const double* updated_exp_ubeta_vals = alt_exp_ubeta_vals + (u_categ * block_n_days);

    for (int r = 0; r < block_n_days; ++r) {

	const double curr_updated_exp_ubeta = updated_exp_ubeta_vals[r];

	block_ubeta_vals[r] = log(curr_updated_exp_ubeta);
	block_exp_ubeta_vals[r] = curr_updated_exp_ubeta;
    }
}




void UGen::update_utau(UProdTau& utau,
		       const int u_categ,
		       const double* alt_utau_vals,
		       const UMissBlock* const miss_block) {

    const int block_beg_sex_idx = miss_block->beg_sex_idx;
    const int block_n_sex_days  = miss_block->n_sex_days;

    double* block_utau_vals         = utau.vals() + block_beg_sex_idx;
    const double* updated_utau_vals = alt_utau_vals + (u_categ * block_n_sex_days);

    std::copy(updated_utau_vals,
	      updated_utau_vals + block_n_sex_days,
	      block_utau_vals);
}
