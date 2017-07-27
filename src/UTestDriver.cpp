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


bool g_record_status;

int g_n_samp;
double g_eps;
Rcpp::NumericVector g_phi_specs;
Rcpp::NumericVector g_test_data_phi;
Rcpp::NumericVector g_test_data_phi_samples;
Rcpp::NumericVector g_xi_vals;
Rcpp::List g_subj_day_blocks;

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
	       int n_samp,
	       Rcpp::NumericVector xi_vals,
	       Rcpp::NumericVector test_data_phi,
	       Rcpp::NumericVector test_data_phi_samples) {

    g_ut_factory = UTestFactory(U,
				X_rcpp,
				w_day_blocks,
				w_to_days_idx,
				w_cyc_to_cyc_idx,
				subj_day_blocks,
				day_to_subj_idx,
				gamma_specs,
				phi_specs,
				fw_len,
				n_burn,
				n_samp,
				xi_vals,
				test_data_phi,
				test_data_phi_samples);


    g_record_status = true;
    g_eps = 0.00000001;

    g_n_samp = n_samp;
    g_phi_specs = phi_specs;
    g_xi_vals = xi_vals;
    g_test_data_phi = test_data_phi;
    g_test_data_phi_samples = test_data_phi_samples;
    g_subj_day_blocks = subj_day_blocks;

    CppUnit::TextUi::TestRunner runner;
    runner.addTest(PhiGenTest::suite());
    runner.addTest(XiGenTest::suite());


    bool out = runner.run();

    return out ? 0 : 1;
}
