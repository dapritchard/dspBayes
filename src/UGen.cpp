#include <cmath>
#include "Rcpp.h"

#include "CoefGen.h"
#include "UGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XiGen.h"

#define REF_CELL     -1
#define IN_NON_PREG_CYC -1
#define NONMISS_SEX_DAY -1




UGen::UGen() :
    m_col_start(0),
    m_col_end(0),
    m_ref_col(0),
    // m_n_categs(0),
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





// void calc_posterior_x_probs() {


// }



// void UGen::calc_posterior_w_probs(double* posterior_w_probs,
// 				  const WGen& W,
// 				  const XiGen& xi,
// 				  const CoefGen coefs,
// 				  const UProdBeta ubeta,
// 				  const UMissBlock* const miss_block) const {

//     // each element of `categ_exp_ubeta_arrs` is used to point to values of
//     // `exp(U * beta)` affected by the missing covariate and corresponding to a
//     // particular choice of covariate category
//     const double* categ_exp_ubeta_arrs[m_n_categs];

//     // storage used for values of `exp(U * beta)` affected by the missing
//     // covariate and for various choices of covariate categories.  We need one
//     // less than the number of categories b/c we already know the values for the
//     // current chioce of category.
//     double alt_exp_ubeta_vals[ miss_block->n_days * (m_n_categs - 1) ];

//     // calculate p(W | U) for the various choices of missing covariate and store
//     // in `posterior_w_probs`
//     calc_ubeta_possibs(alt_exp_ubeta_vals, categ_exp_ubeta_arrs, coefs, ubeta, miss_block);
//     calc_posterior_w_from_ubetas(posterior_w_probs, W, xi, miss_block, categ_exp_ubeta_arrs);
// }




void UGen::calc_ubeta_possibs(double* alt_exp_ubeta_vals,
			      const double** categ_exp_ubeta_arrs,
			      const CoefGen& coefs,
			      const UProdBeta& ubeta,
			      const UMissBlock* const miss_block) const {

    const int block_beg_day_idx = miss_block->beg_day_idx;
    const int block_n_days      = miss_block->n_days;
    const int block_u_col       = miss_block->u_col;

    const double* block_exp_ubeta_vals = ubeta.exp_vals()[block_beg_day_idx];
    const double* block_w_idx = m_vals()[block_w_day_idx];

    const double* beta_coefs = coefs.vals();
    const int* w_vals        = W.vals();
    const double xi_i        = xi.vals()[miss_block->subj_idx];

    //     const int* w_vals      = W.vals();

    // each iteration calculates `p(W | U, data)` for the current category of
    // the missing covariate
    for (int j = m_col_start; j < m_col_end; ++j) {

	double exp_beta_diff;
	double log_dpois_sum;

	// the value of `exp(beta_j - beta^{*})` where `beta^{*}` is the
	// coefficient of beta that corresponds to the category of `U` chosen
	// for the previous sample.  Either `beta_j` or `beta^{*}` or both may
	// correspond the the reference cell category, in which case the
	// corresponding value is implicitely 0.

	if (j == block_u_col) {
	    exp_beta_diff = 1;
	}
	else {
	    const double beta_j = (j == m_ref_col) ? 0.0 : beta_coefs[j];
	    exp_beta_diff = (block_u_col == REF_CELL) ?
		exp(beta_j):
		exp(beta_j - beta_coefs[block_u_col]);
	}

	for (int r = 0 ; r < block_n_days; ++r) {

		// the mean value for the Poisson distribution of `W_ijk`, and the
		// index in `W.vals()` corresponding to observation ijk
		const double curr_mean_val = xi_i * block_exp_ubeta_vals[r] * exp_beta_diff;
		const int curr_w_idx       = block_w_idx[r];

		// the value of `p(W_ijk | U, data)`.  If `curr_w_idx` has the
		// flag value indicating that observation `ijk` corresponds to a
		// day in a non-pregnancy cycle, then `W_ijk` has a value of 0.
		log_dpois_sum += (curr_w_idx == IN_NON_PREG_CYC) ?
		    -curr_mean_val :
		    R::dpois(w_vals[curr_w_idx], mean_val, 1);
	}

	// exponentiate and store the result
	*posterior_w_probs++ = exp(log_dpois_sum);

	// // case: the j-th column is the one for which the current `U` has a
	// // value of 1.  No need to calculate `exp(U * beta)` for this choice of
	// // `U`.
	// if (j == block_u_col) {
	//     *categ_exp_ubeta_arrs++ = ubeta_exp_vals + block_beg_idx;
	// }

	// // case: the j-th column is not the one for which the current `U` has a
	// // value of 1.  Need to calculate the value `exp(U * beta)` for the
	// // modified choice of `U`.
	// else {

	//     // the value of `exp(beta_j - beta^{*})` where `beta^{*}` is the
	//     // coefficient of beta that corresponds to the category of `U`
	//     // chosen for the previous sample.  Either `beta_j` or `beta^{*}` or
	//     // both may correspond the the reference cell category, in which
	//     // case the corresponding value is implicitely 0.
	//     const double beta_j = (j == m_ref_col) ?
	// 	0.0 :
	// 	beta_coefs[j];
	//     const double exp_beta_diff = (block_u_col == REF_CELL) ?
	// 	exp(beta_j):
	// 	exp(beta_j - beta_coefs[block_u_col]);
	//     // const double exp_beta_diff = (curr_u_col == REF_CELL) ?
	//     // 	exp(coefs_vals[j]) :
	//     // 	exp(coefs_vals[j] - coefs_vals[curr_u_col]);

	//     // step through the observations affected by the missing covariate
	//     // and calculate `exp(U * beta)` for the modified choice of `U`
	//     for (int r = 0 ; r < curr_n_days; ++r) {
	// 	alt_exp_ubeta_vals[r] = ubeta_exp_vals[curr_beg_idx + r] * exp_beta_diff;
	//     }

	//     // point to the modified values of `exp(U * beta)` corresponding to
	//     // a 1 for the current column
	//     *categ_exp_ubeta_arrs++ = alt_exp_ubeta_vals;
	//     // point to the next free location in the storage for the modified
	//     // values
	//     alt_exp_ubeta_vals += curr_n_days;
	// }
    }

    // // case: the categorical variable corresponding to the current group of
    // // columns uses reference-cell coding, so we have to obtain the values of
    // // `exp(U * beta)` corresponding to the reference cell
    // if (m_is_ref_cell_coding) {

    // 	// case: the currrent value of `U` corresponds to the reference cell
    // 	if (curr_u_col == REF_CELL) {
    // 	    *categ_exp_ubeta_arrs = ubeta_exp_vals + curr_beg_idx;
    // 	}

    // 	// case: the reference cell is not the one for which the previous value
    // 	// of `U` has a value of 1.  Need to calculate the value `exp(U * beta)`
    // 	// for the modified choice of `U`.
    // 	else {

    // 	    // the exponentiated difference between the beta coefficient
    // 	    // corresponding to the reference cell and the beta coefficient
    // 	    // corresponding to column that currently has a value of 1.
    // 	    const double exp_beta_diff = std::exp(- coefs_vals[curr_u_col]);

    // 	    // step through the observations affected by the missing covariate
    // 	    // and calclate `U * beta` for the modified choice of `U`
    // 	    for (int r = 0 ; r < curr_n_days; ++r) {
    // 		alt_exp_ubeta_vals[r] = ubeta_exp_vals[curr_beg_idx + r] * exp_beta_diff;
    // 	    }

    // 	    // point to the modified values of `exp(U * beta)` corresponding to
    // 	    // a 1 for the current column
    // 	    *categ_exp_ubeta_arrs = alt_exp_ubeta_vals;
    // 	}
    // }
}




// void UGen::calc_posterior_w_from_utaus(double* posterior_w_probs,
// 				       const WGen& W,
// 				       const XiGen& xi,
// 				       const UMissBlock* const miss_block,
// 				       const double* const * categ_exp_ubeta_arrs) const {

//     const int curr_beg_idx = miss_block->beg_idx;
//     const int curr_n_days  = miss_block->n_days;
//     const double xi_i      = xi.vals()[miss_block->subj_idx];
//     const int* w_vals      = W.vals();

//     // each iteration calculates `p(W | U, data)` for the current missing
//     // covariate
//     for (int j = 0; j < m_n_categs; ++j) {

// 	// `curr_exp_ubeta_vals` points to the values of `exp(U * beta)` for the
// 	// choice of U corresponding to the j-th category
// 	const double* curr_exp_ubeta_vals = categ_exp_ubeta_arrs[j];
// 	double log_dpois_sum = 0.0;

// 	// each iteration adds the value of `log p(W_ijk | U, data)` for the for
// 	// the choice of U corresponding to the j-th category to
// 	// `log_dpois_sum`.
// 	for (int r = 0; r < curr_n_days; ++r) {

// 	    // the mean value for the Poisson distribution of `W_ijk`, and the
// 	    // index in `W.vals()` corresponding to observation ijk
// 	    double mean_val = xi_i * curr_exp_ubeta_vals[r];
// 	    int curr_w_idx  = m_w_idx[curr_beg_idx + r];

// 	    // the value of `p(W_ijk | U, data)`.  If `curr_w_idx` has the flag
// 	    // value indicating that observation ijk corresponds to a day in a
// 	    // non-pregnancy cycle, then `W_ijk` has a value of 0.
// 	    log_dpois_sum += (curr_w_idx == IN_NON_PREG_CYC) ?
// 		-mean_val :
// 		R::dpois(w_vals[curr_w_idx], mean_val, 1);
// 	}

// 	// exponentiate and store the result
// 	*posterior_w_probs++ = exp(log_dpois_sum);
//     }
// }




// void UGen::calc_alt_utaus(double* alt_utau_vals,
// 			  const UProdBeta& utau,
// 			  const UMissBlock* const miss_block) const {

//     const int curr_beg_idx    = miss_block->beg_idx;
//     const int curr_n_sex_days = miss_block->n_sex_days;
//     const int curr_u_col      = miss_block->u_col;

//     const XGen::XMissDay* x_miss_day = X.miss_day();
//     const double sex_coef            = X.sex_coef();

//     const double* utau_coefs = utau.coefs();
//     const double* utau_vals  = utau.vals();

//     // each iteration obtains `U * tau` corresponding to a value of 1 for the
//     // dummy coding in the j-th variable for the observations that are affected
//     // by the missing covariate
//     for (int j = m_col_start; j < m_col_end; ++j) {

// 	// the difference between tau_j and the tau value corresponding to the
// 	// current value of `U`.  The current value is implicitely 0 if the
// 	// current category for the missing covariate corresponds to the
// 	// reference cell.
// 	const double tau_diff = (curr_u_col == REF_CELL) ?
// 	    utau_coefs[j] :
// 	    utau_coefs[j] - utau_coefs[curr_u_col];

// 	// step through the observations affected by the missing covariate and
// 	// calculate `U * tau` for the modified choice of `U`
// 	for (int r = 0; r < curr_n_days; ++r) {

// 	    // only calculate `U * tau` for days in the intercourse variable was
// 	    // missing, i.e. days in which `P(X | U)` is a random variable
// 	    const int curr_x_day_idx = m_miss_day[curr_beg_sex_idx + r];
// 	    if (curr_x_day_idx == NONMISS_SEX_DAY) {
// 		continue;
// 	    }

// 	    // a value of <= 0 corresponds to no sex on the previous day, and >=
// 	    // 1 corresponds to yes sex on the previous day
// 	    *alt_utau_vals++ = (x_miss_day[curr_x_day_idx] <= 0) ?
// 		utau_vals[curr_x_day_idx] :
// 		utau_vals[curr_x_day_idx] + sex_coef;
// 	}
//     }

//     // case: the categorical variable corresponding to the current group of
//     // columns uses reference-cell coding, so we have to obtain the values of `U
//     // * tau` corresponding to the reference cell
//     if (m_is_ref_cell_coding) {

// 	// the difference between the tau coefficient corresponding to the
// 	// reference cell and the tau coefficient corresponding to column
// 	// that currently has a value of 1
// 	const double tau_diff = (curr_u_col == REF_CELL) ?
// 	    0.0 :
// 	    -utau_coefs[curr_u_col];

// 	// step through the observations affected by the missing covariate
// 	// and calculate `U * tau` for the modified choice of `U`
// 	for (int r = 0; r < curr_n_days; ++r) {

// 	    // only calculate `U * tau` for days in the intercourse variable was
// 	    // missing, i.e. days in which `P(X | U)` is a random variable
// 	    const int curr_x_day_idx = m_miss_day[curr_beg_idx + r];
// 	    if (curr_x_day_idx == NONMISS_SEX_DAY) {
// 		continue;
// 	    }

// 	    // a value of <= 0 corresponds to no sex on the previous day, and >=
// 	    // 1 corresponds to yes sex on the previous day
// 	    *alt_utau_vals++ = (x_miss_day[curr_x_day_idx].prev <= 0) ?
// 		utau_vals[curr_x_day_idx] :
// 		utau_vals[curr_x_day_idx] + sex_coef;
// 	}
//     }
// }




void UGen::calc_posterior_x(double* posterior_x_probs,
			    const XGen& X,
			    const UProdTau& utau,
			    const UMissBlock* const miss_block) const {

    const int block_beg_sex_idx = miss_block->beg_sex_idx;
    const int block_n_sex_days  = miss_block->n_sex_days;
    const int block_u_col       = miss_block->u_col;

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

	// step through the observations affected by the missing covariate and
	// calculate `U * tau` for the modified choice of `U`
	for (int r = 0; r < block_n_sex_days; ++r) {

	    const int curr_x_day_idx = m_x_idx[block_beg_sex_idx + r];
	    const int curr_day_idx   = x_miss_day[curr_x_day_idx].idx;
	    const int curr_sex_prev  = x_miss_day[curr_x_day_idx].prev;

	    // // only calculate `U * tau` for days in the intercourse variable was
	    // // missing, i.e. days in which `P(X | U)` is a random variable
	    // const int curr_x_day_idx = m_miss_day[block_beg_idx + r].idx;
	    // if (curr_x_day_idx == NONMISS_SEX_DAY) {
	    // 	continue;
	    // }

	    // a value of <= 0 corresponds to no sex on the previous day, and >=
	    // 1 corresponds to yes sex on the previous day
	    double curr_utau_val = (curr_sex_prev <= 0) ?
		utau_vals[curr_x_day_idx] + tau_diff :
		utau_vals[curr_x_day_idx] + tau_diff + sex_coef;

	    // add `-log P(X_ijk | U_ijk)` to the running total, where `P(X_ijk
	    // | U_ijk)` follows a logistic regression model
	    neg_log_probs_sum += (curr_day_idx) ?
		log(1 + exp(- curr_utau_val)) :
		log(1 + exp(curr_utau_val));
	}

	// `P(X = x | U) = exp( sum_{ijk: ijk is affected by U} log( P(X_ijk |
	// U_ijk) ) )`
	*posterior_x_probs++ = exp(- neg_log_probs_sum);
    }

    // // case: the categorical variable corresponding to the current group of
    // // columns uses reference-cell coding, so we have to obtain the values of `U
    // // * tau` corresponding to the reference cell
    // if (m_is_ref_cell_coding) {

    // 	// the difference between the tau coefficient corresponding to the
    // 	// reference cell and the tau coefficient corresponding to column
    // 	// that currently has a value of 1
    // 	const double tau_diff = (curr_u_col == REF_CELL) ?
    // 	    0.0 :
    // 	    -utau_coefs[curr_u_col];

    // 	// step through the observations affected by the missing covariate
    // 	// and calculate `U * tau` for the modified choice of `U`
    // 	for (int r = 0; r < curr_n_days; ++r) {

    // 	    // only calculate `U * tau` for days in the intercourse variable was
    // 	    // missing, i.e. days in which `P(X | U)` is a random variable
    // 	    const int curr_x_day_idx = m_miss_day[curr_beg_idx + r];
    // 	    if (curr_x_day_idx == NONMISS_SEX_DAY) {
    // 		continue;
    // 	    }

    // 	    // a value of <= 0 corresponds to no sex on the previous day, and >=
    // 	    // 1 corresponds to yes sex on the previous day
    // 	    *alt_utau_vals++ = (x_miss_day[curr_x_day_idx].prev <= 0) ?
    // 		utau_vals[curr_x_day_idx] :
    // 		utau_vals[curr_x_day_idx] + sex_coef;

    // 	    // add `-log P(X_ijk | U_ijk)` to the running total
    // 	    neg_log_probs_sum += ( x_vals[ x_miss_day[curr_x_day_idx].idx ] ) ?
    // 		log(1 + exp(- curr_utau_val)) :
    // 		log(1 + exp(curr_utau_val));
    // 	}
    // }
}




// void UGen::calc_posterior_x_from_utaus(double* posterior_w_probs,
// 				       const XiGen& xi,
// 				       const XGen& X,
// 				       const UMissBlock* const miss_block,
// 				       const double* alt_utau_vals) const {

//     const int curr_beg_idx = miss_block->beg_idx;
//     const int curr_n_days  = miss_block->n_days;
//     const double xi_i      = xi.vals()[miss_block->subj_idx];
//     const int* x_vals      = X.vals();

//     // each iteration calculates `p(X | U)` for the current choice of `U`
//     for (int j = 0; j < m_n_categs; ++j) {

// 	double neg_log_probs_sum = 0.0;

// 	// each iteration adds the value of `-log p(X_ijk | U, data)` for the
// 	// for the choice of U corresponding to the j-th category to
// 	// `neg_log_probs_sum`.
// 	for (int r = 0; r < curr_n_days; ++r) {

// 	    // // the mean value for the Poisson distribution of `W_ijk`, and the
// 	    // // index in `W.vals()` corresponding to observation ijk
// 	    // double mean_val = xi_i * curr_exp_ubeta_vals[r];
// 	    // int curr_w_idx  = m_w_idx[curr_beg_idx + r];

// 	    // only calculate `P(X | U)` for days in the intercourse variable
// 	    // was missing, i.e. days in which `P(X | U)` is a random variable
// 	    const int curr_x_day_idx = m_x_idx[curr_beg_idx + r];
// 	    if (curr_x_day_idx == NONMISS_SEX_DAY) {
// 		continue;
// 	    }

// 	    neg_log_probs_sum += ( x_vals[ x_miss_day[curr_x_day_idx].idx ] ) ?
// 		log(1 + exp(- curr_utau_val)) :
// 		log(1 + exp(curr_utau_val));

// 	    // // the value of `p(W_ijk | U, data)`.  If `curr_w_idx` has the flag
// 	    // // value indicating that observation ijk corresponds to a day in a
// 	    // // non-pregnancy cycle, then `W_ijk` has a value of 0.
// 	    // log_dpois_sum += (curr_w_idx == IN_NON_PREG_CYC) ?
// 	    // 	-mean_val :
// 	    // 	R::dpois(w_vals[curr_w_idx], mean_val, 1);
// 	}

// 	// exponentiate and store the result
// 	*posterior_w_probs = std::exp(log_dpois_sum);
//     }
// }
