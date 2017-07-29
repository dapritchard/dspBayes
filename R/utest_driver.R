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

    # sample women-specific fecundability multiplier
    set.seed(100L)
    xi <- rgamma(length(dsp_data$subj_day_blocks), 1, 1)

    # sample latent pregnancy variable W
    set.seed(101L)
    W <- rpois(length(dsp_data$w_to_days_idx), 1.5)

    # sample U * beta
    set.seed(102L)
    ubeta <- rnorm(length(dsp_data$X))

    # phi testing data
    phi_seed <- 99L
    out_utest_phi <- utest_phi(xi, seed_val = phi_seed)

    # W testing data
    w_seed <- 207L
    target_samples_w <- utest_w(dsp_data, xi, ubeta, w_seed)

    # xi testing data
    xi_seed <- 21L
    target_samples_xi <- utest_xi(dsp_data, W, ubeta, phi_specs["mean"], xi_seed)

    # collect seeds
    seed_vals <- c(phi = phi_seed,
                   w   = w_seed,
                   xi  = xi_seed)

    # collect testing objects
    test_data <- list(input_ubeta        = ubeta,
                      input_w            = W,
                      input_xi           = xi,
                      target_data_phi    = out_utest_phi$target_data,
                      target_samples_phi = out_utest_phi$target_samples,
                      target_samples_w   = target_samples_w,
                      target_samples_xi  = target_samples_xi,
                      seed_vals          = seed_vals,
                      epsilon            = 1e-12)

    n_samp <- 100    #  ************  TODO remove this line  ***********************

    # pass testing data to C++ testing driver
    utest_cpp_(U                 = dsp_data$U,
               X_rcpp            = dsp_data$X,
               w_day_blocks      = dsp_data$w_day_blocks,
               w_to_days_idx     = dsp_data$w_to_days_idx,
               w_cyc_to_subj_idx = dsp_data$w_cyc_to_subj_idx,
               subj_day_blocks   = dsp_data$subj_day_blocks,
               day_to_subj_idx   = dsp_data$day_to_subj_idx,
               gamma_specs       = gamma_hyper_list,
               phi_specs         = phi_specs,
               fw_len            = 5,
               n_burn            = 0,
               n_samp            = n_samp,
               test_data         = test_data)
}
