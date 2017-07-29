#include <algorithm>
#include "Rcpp.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestWGen.h"

using Rcpp::as;

extern UTestFactory g_ut_factory;




WGenTest::WGenTest() :
    target_samples(g_ut_factory.target_samples_w),
    n_preg_days(g_ut_factory.w_to_days_idx.size()),
    n_preg_cyc(g_ut_factory.preg_cyc.size()),
    seed_val(as<double>(g_ut_factory.seed_vals["W"])),
    epsilon(UTestFactory::epsilon) {
}




void WGenTest::setUp() {

    // construct W
    W = g_ut_factory.W();

    // construct xi
    xi = g_ut_factory.xi_no_rec();

    // constuct ubeta
    ubeta = g_ut_factory.ubeta();
}




void WGenTest::tearDown() {
    delete W;
    delete xi;
    delete ubeta;
}




void WGenTest::test_constructor() {
    CPPUNIT_ASSERT_EQUAL(n_preg_days, W->m_n_preg_days);
    CPPUNIT_ASSERT_EQUAL(n_preg_cyc, W->m_n_preg_cyc);
}




void WGenTest::test_sample() {

    W->sample(*xi, *ubeta);

    // check values of samples
    CPPUNIT_ASSERT(std::equal(target_samples.begin(),
			      target_samples.end(),
			      W->m_vals,
			      UTestFactory::eq_dbl));

    // check values of `sum_jk W_ijk`
    int w_ctr = 0;
    int* w_vals = W->m_vals;
    int* w_sums = W->m_sums;
    for (const PregCyc* curr = W->m_preg_cyc; curr < W->m_preg_cyc + W->m_n_preg_cyc; ++curr) {
    	int sum_val = std::accumulate(w_vals + w_ctr, w_vals + w_ctr + curr->n_days, 0.0);
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(sum_val, *w_sums++, epsilon);
    	w_ctr += curr->n_days;
    }
}
