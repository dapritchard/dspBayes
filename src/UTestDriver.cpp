#include "Rcpp.h"
#include "CoefGen.h"
#include "PhiGen.h"
#include "WGen.h"
#include "XiGen.h"

#include "cppunit/ui/text/TestRunner.h"
#include "UTestFactory.h"
#include "UTestPhiGen.h"
#include "UTestXiGen.h"

extern int* d2s;

UTestFactory g_ut_factory;

// preg_cycle            used when sampling W
// w_to_days_idx         categorical gamma: a_tilde
// w_cyc_to_subj_idx     used when sampling xi (first term)
// fw_len                how much memory to set aside when sampling W in a cycle
// subj_day_block        used when sampling xi (second term)
// gamma_specs           gamma hyperparameters
// phi_specs             phi hyperparameters




// [[Rcpp::export]]
int utest_cpp_(Rcpp::NumericMatrix U,
	       Rcpp::IntegerVector X_rcpp,
	       Rcpp::List w_day_blocks,
	       Rcpp::IntegerVector w_to_days_idx,
	       Rcpp::IntegerVector w_cyc_to_subj_idx,
	       Rcpp::List subj_day_blocks,
	       Rcpp::IntegerVector day_to_subj_idx,
	       Rcpp::List gamma_specs,
	       Rcpp::NumericVector phi_specs,
	       int fw_len,
	       int n_burn,
	       int n_samp,
	       Rcpp::List test_data) {

    UTestFactory::epsilon = Rcpp::as<double>(test_data["epsilon"]);
    g_ut_factory = UTestFactory(U,
				X_rcpp,
				w_day_blocks,
				w_to_days_idx,
				w_cyc_to_subj_idx,
				subj_day_blocks,
				day_to_subj_idx,
				gamma_specs,
				phi_specs,
				fw_len,
				n_burn,
				n_samp,
				test_data);

    CppUnit::TextUi::TestRunner runner;
    runner.addTest(PhiGenTest::suite());
    runner.addTest(XiGenTest::suite());

    // run() returns true if successful, false otherwise
    return runner.run() ? 0 : 1;
}
