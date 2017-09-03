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
    u_prior_probs(as<Rcpp::NumericVector>(as<Rcpp::List>(g_ut_factory.u_miss_info[var_idx])["u_prior_probs"])),
    input_w_probs(as<Rcpp::NumericVector>(g_ut_factory.input_u_categ["posterior_w_probs"])),
    input_x_probs(as<Rcpp::NumericVector>(g_ut_factory.input_u_categ["posterior_x_probs"])),
    input_block_idx(g_ut_factory.target_data_u_categ["block_idx"]),
    target_x_probs(as<Rcpp::IntegerVector>(g_ut_factory.target_samples_u_categ["posterior_x"])),
    target_sample_covs(as<Rcpp::IntegerVector>(g_ut_factory.target_samples_u_categ["sample_covs"])),
    target_alt_utau_vals(as<Rcpp::NumericVector>(g_ut_factory.target_samples_u_categ["alt_utau_vals"])),
    seed_val(g_ut_factory.seed_vals["u_categ"]),
    epsilon(UTestFactory::epsilon)
{}




void UGenVarCategTest::setUp() {

    u_var = g_ut_factory.u_categ();




    // W     = g_ut_factory.W();
    // xi    = g_ut_factory.xi_no_rec();
    // ubeta = g_ut_factory.ubeta();
    utau  = g_ut_factory.utau();
    X     = g_ut_factory.X();

    // // copy X data so that tests don't cause persistent changes
    // x_rcpp_copy = new Rcpp::IntegerVector(X_rcpp.begin(), X_rcpp.end());
    // X = new XGen(*x_rcpp_copy, miss_cyc_rcpp, miss_day_rcpp, cohort_sex_prob, sex_coef);

}




void UGenVarCategTest::tearDown() {
    delete u_var;
    delete X;
    // delete W;
    // delete xi;
    // delete ubeta;
    delete utau;
    // delete x_rcpp_copy;
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
    CPPUNIT_ASSERT(std::equal(u_prior_probs.begin(), u_prior_probs.end(), u_var->m_u_prior_probs));
    // TODO: test `m_miss_block`?
}




void UGenVarCategTest::test_calc_posterior_x() {

    double posterior_x_probs[u_var->m_n_categs];
    double alt_utau_vals[u_var->m_max_n_sex_days_miss];

    u_var->calc_posterior_x(posterior_x_probs,
    			    alt_utau_vals,
    			    *X,
    			    *utau,
    			    u_var->m_miss_block + input_block_idx);

    CPPUNIT_ASSERT(std::equal(target_x_probs.begin(),
    			      target_x_probs.end(),
    			      posterior_x_probs));

    CPPUNIT_ASSERT(std::equal(target_alt_utau_vals.begin(),
    			      target_alt_utau_vals.end(),
    			      alt_utau_vals));
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
