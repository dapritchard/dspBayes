utest_cpp <- function(dsp_data,
                      nSamp      = 1e4,
                      nBurn      = 5e3,
                      nThin      = 1L,
                      hypGam     = NULL,
                      tuningGam  = NULL,
                      hypPhi     = NULL,
                      tuningPhi  = 0.3,
                      trackProg  = "percent",
                      progQuants = seq(0.1, 1.0, 0.1)) {

    # stub functions for gamma and phi specs
    gamma_hyper_list <- get_gamma_specs(dsp_data)
    phi_specs <- get_phi_specs()

    n_samp <- 100

    out <- utest_cpp_(U                = dsp_data$U,
                       X_rcpp           = dsp_data$X,
                       w_day_blocks     = dsp_data$w_day_blocks,
                       w_to_days_idx    = dsp_data$w_to_days_idx,
                       w_cyc_to_cyc_idx = dsp_data$w_cyc_to_cyc_idx,
                       subj_day_blocks  = dsp_data$subj_day_blocks,
                       day_to_subj_idx  = dsp_data$day_to_subj_idx,
                       gamma_specs      = gamma_hyper_list,
                       phi_specs        = phi_specs,
                       fw_len           = 5,
                       n_burn           = 0,
                       n_samp           = n_samp)


}
