#include "Rcpp.h"
#include "CoefGen.h"
#include "PhiGen.h"
#include "WGen.h"
#include "XiGen.h"

#include "cppunit/ui/text/TestRunner.h"
#include "UTestPhiGen.h"

extern int* d2s;

Rcpp::NumericVector g_phi_specs;
int g_n_samp;

// preg_cycle            used when sampling W
// w_to_days_idx         categorical gamma: a_tilde
// w_cyc_to_cyc_idx      used when sampling xi (first term)
// fw_len                how much memory to set aside when sampling W in a cycle
// subj_day_block        used when sampling xi (second term)
// gamma_specs           gamma hyperparameters
// phi_specs             phi hyperparameters




// [[Rcpp::export]]
int utest_cpp_(Rcpp::NumericMatrix U,
	       Rcpp::IntegerVector X_rcpp,
	       Rcpp::List w_day_blocks,
	       Rcpp::IntegerVector w_to_days_idx,
	       Rcpp::IntegerVector w_cyc_to_cyc_idx,
	       Rcpp::List subj_day_blocks,
	       Rcpp::IntegerVector day_to_subj_idx,
	       Rcpp::List gamma_specs,
	       Rcpp::NumericVector phi_specs,
	       int fw_len,
	       int n_burn,
	       int n_samp) {

    g_phi_specs = phi_specs;
    g_n_samp = n_samp;

    CppUnit::TextUi::TestRunner runner;
    runner.addTest(PhiGenTest::suite());

    bool out = runner.run();

    return out ? 0 : 1;
}
