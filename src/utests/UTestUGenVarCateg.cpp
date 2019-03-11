#include <algorithm>
#include "Rcpp.h"

#include "CoefGen.h"
#include "UGenVar.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "UTestFactory.h"
#include "UTestUGenVarCateg.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"


// using Rcpp::IntegerVector;
// using Rcpp::NumericVector;
using Rcpp::as;

extern UTestFactory g_ut_factory;


UGenVarCategTest::UGenVarCategTest() :
    var_idx(g_ut_factory.target_data_u_categ["var_idx"]),
    u_rcpp(g_ut_factory.u_rcpp),
    n_days(g_ut_factory.target_data_u_categ["n_days"]),
    w_idx(g_ut_factory.u_preg_map),
    x_idx(g_ut_factory.u_sex_map),
    col_start(g_ut_factory.target_data_u_categ["col_start"]),
    col_end(g_ut_factory.target_data_u_categ["col_end"]),
    ref_col(g_ut_factory.target_data_u_categ["ref_col"]),
    n_categs(g_ut_factory.target_data_u_categ["n_categs"]),
    max_n_days_miss(g_ut_factory.target_data_u_categ["max_n_days_miss"]),
    max_n_sex_days_miss(g_ut_factory.target_data_u_categ["max_n_sex_days_miss"]),
    log_u_prior_probs(as<Rcpp::NumericVector>(as<Rcpp::List>(g_ut_factory.u_miss_info[var_idx])["log_u_prior_probs"])),
    input_w_probs(as<Rcpp::NumericVector>(g_ut_factory.input_u_categ["posterior_w_probs"])),
    input_x_probs(as<Rcpp::NumericVector>(g_ut_factory.input_u_categ["posterior_x_probs"])),
    input_block_idx(g_ut_factory.target_data_u_categ["block_idx"]),
    target_w_probs(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["posterior_w"])),
    target_x_probs(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["posterior_x"])),
    target_sample_covs(as<Rcpp::IntegerVector>(g_ut_factory.target_samples_u_categ["sample_covs"])),
    target_alt_exp_ubeta_vals(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["alt_exp_ubeta_vals"])),
    target_alt_utau_vals(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["alt_utau_vals"])),
    target_categ_update(as<Rcpp::IntegerVector>(g_ut_factory.target_samples_u_categ["categ_update"])),
    target_ubeta_update(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["ubeta_update"])),
    target_utau_update(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["utau_update"])),
    target_u_update(as<Rcpp::NumericMatrix>(g_ut_factory.target_samples_u_categ["u_update"])),
    seed_val(g_ut_factory.seed_vals["u_categ"]),
    epsilon(UTestFactory::epsilon)
{}




void UGenVarCategTest::setUp() {

    W     = g_ut_factory.W();
    xi    = g_ut_factory.xi_no_rec();
    coefs = g_ut_factory.coefs();
    ubeta = g_ut_factory.ubeta();
    utau  = g_ut_factory.utau();
    X     = g_ut_factory.X();

    // so changes to underlying data in `u_var` aren't persistent
    u_rcpp_copy = new Rcpp::NumericMatrix;
    *u_rcpp_copy = clone(u_rcpp);
    u_var = g_ut_factory.u_categ(u_rcpp_copy);

    // so changes to underlying data in `utau` aren't persistent
    utau_vals_copy = new double[utau->m_vals_rcpp.size()];
    std::copy(utau->m_vals_rcpp.begin(), utau->m_vals_rcpp.end(), utau_vals_copy);
    utau->m_vals = utau_vals_copy;
}




void UGenVarCategTest::tearDown() {

    delete u_var;
    delete W;
    delete xi;
    delete coefs;
    delete ubeta;
    delete utau;
    delete X;

    delete u_rcpp_copy;
    delete[] utau_vals_copy;
}




void UGenVarCategTest::test_constructor() {

    // parent class
    CPPUNIT_ASSERT(std::equal(u_rcpp.begin() + ((int) u_rcpp.nrow() * col_start),
                              u_rcpp.begin() + ((int) u_rcpp.nrow() * ref_col),
                              u_var->m_u_var_col));
    CPPUNIT_ASSERT_EQUAL(n_days, u_var->m_n_days);
    CPPUNIT_ASSERT(std::equal(w_idx.begin(), w_idx.end(), u_var->m_w_idx));
    CPPUNIT_ASSERT(std::equal(x_idx.begin(), x_idx.end(), u_var->m_x_idx));

    // child class
    CPPUNIT_ASSERT_EQUAL(col_start, u_var->m_col_start);
    CPPUNIT_ASSERT_EQUAL(col_end, u_var->m_col_end);
    CPPUNIT_ASSERT_EQUAL(ref_col, u_var->m_ref_col);
    CPPUNIT_ASSERT_EQUAL(n_categs, u_var->m_n_categs);
    CPPUNIT_ASSERT_EQUAL(max_n_days_miss, u_var->m_max_n_days_miss);
    CPPUNIT_ASSERT_EQUAL(max_n_sex_days_miss, u_var->m_max_n_sex_days_miss);
    CPPUNIT_ASSERT(std::equal(log_u_prior_probs.begin(), log_u_prior_probs.end(), u_var->m_log_u_prior_probs));
    // TODO: test `m_miss_block`?
}




void UGenVarCategTest::test_sample() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // sample missing covariate values
    u_var->sample(*W, *xi, *coefs, *X, *ubeta, *utau);

    // check sampled covariate
    for (int i = 0; i < target_categ_update.size(); ++i) {
        CPPUNIT_ASSERT_EQUAL(target_categ_update[i], u_var->m_miss_block[i].u_col);
    }

    // check `U * beta` update
    CPPUNIT_ASSERT(std::equal(target_ubeta_update.begin(),
                              target_ubeta_update.end(),
                              ubeta->vals(),
                              UTestFactory::eq_dbl));

    // TODO: test exp(ubeta) ?

    // check `U * tau` update
    CPPUNIT_ASSERT(std::equal(target_utau_update.begin(),
                              target_utau_update.end(),
                              utau->vals(),
                              UTestFactory::eq_dbl));

    // check `U` update
    CPPUNIT_ASSERT(std::equal(target_u_update.begin(),
                              target_u_update.end(),
                              u_var->m_u_var_col,
                              UTestFactory::eq_dbl));
}




void UGenVarCategTest::test_calc_log_condit_w() {

    double log_condit_w_probs[u_var->m_n_categs];
    double alt_exp_ubeta_vals[u_var->m_max_n_days_miss * u_var->m_n_categs];

    u_var->calc_log_condit_w(log_condit_w_probs,
                             alt_exp_ubeta_vals,
                             *W,
                             *xi,
                             *coefs,
                             *X,
                             *ubeta,
                             u_var->m_miss_block + input_block_idx);

    CPPUNIT_ASSERT(std::equal(target_alt_exp_ubeta_vals.begin(),
                              target_alt_exp_ubeta_vals.end(),
                              alt_exp_ubeta_vals,
                              UTestFactory::eq_dbl));

    CPPUNIT_ASSERT(std::equal(target_w_probs.begin(),
                              target_w_probs.end(),
                              log_condit_w_probs,
                              UTestFactory::eq_dbl));
}




void UGenVarCategTest::test_calc_log_condit_x() {

    double log_condit_x_probs[u_var->m_n_categs];
    double alt_utau_vals[u_var->m_max_n_sex_days_miss * u_var->m_n_categs];

    u_var->calc_log_condit_x(log_condit_x_probs,
                            alt_utau_vals,
                            *X,
                            *utau,
                            u_var->m_miss_block + input_block_idx);

    CPPUNIT_ASSERT(std::equal(target_alt_utau_vals.begin(),
                              target_alt_utau_vals.end(),
                              alt_utau_vals,
                              UTestFactory::eq_dbl));

    CPPUNIT_ASSERT(std::equal(target_x_probs.begin(),
                              target_x_probs.end(),
                              log_condit_x_probs,
                              UTestFactory::eq_dbl));
}




void UGenVarCategTest::test_sample_covariate() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    for (int i = 0; i < target_sample_covs.size(); ++i) {
        CPPUNIT_ASSERT_EQUAL(target_sample_covs[i],
                             u_var->sample_covariate(input_w_probs.begin(),
                                                     input_x_probs.begin()));
    }
}




// void UGenVarCategTest::test_update_ubeta() {

//     // construct target `exp(U * beta)`
//     Rcpp::NumericVector target_exp_ubeta_vals( exp(target_ubeta_vals) );

//     // set the data to the last category and update `ubeta`
//     u_var->update_ubeta(*ubeta,
//                      n_categs - 1,
//                      input_uprod_vals.begin(),
//                      u_var->m_miss_block + input_block_idx);

//     CPPUNIT_ASSERT(std::equal(target_ubeta_vals.begin(),
//                            target_ubeta_vals.end(),
//                            ubeta->m_vals,
//                            UTestFactory::eq_dbl));

//     CPPUNIT_ASSERT(std::equal(target_exp_ubeta_vals.begin(),
//                            target_exp_ubeta_vals.end(),
//                            ubeta->m_exp_vals,
//                            UTestFactory::eq_dbl));
// }
