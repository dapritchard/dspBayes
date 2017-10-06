utest_cpp <- function(dsp_data,
                      n_samp     = 1e4,
                      nBurn      = 5e3,
                      nThin      = 1L,
                      hypGam     = NULL,
                      tuningGam  = NULL,
                      hypPhi     = NULL,
                      tuningPhi  = 0.3,
                      trackProg  = "percent",
                      progQuants = seq(0.1, 1.0, 0.1)) {

    # we've added a sentinal value to the end of `w_to_day_idx`, so we need to
    # remove it before performing the testing routines
    w_to_days_idx_old <- dsp_data$w_to_days_idx
    dsp_data$w_to_days_idx <- w_to_days_idx_old[-length(w_to_days_idx_old)]

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
    ubeta <- rnorm(length(dsp_data$intercourse$X))

    # phi testing data
    phi_seed <- 99L
    out_utest_phi <- utest_phi(xi, seed_val = phi_seed)

    # W testing data
    w_seed <- 207L
    target_samples_w <- utest_w(dsp_data, xi, ubeta, w_seed)

    # X testing data
    x_seed <- 7201L
    out_utest_x <- utest_x(dsp_data, W, xi, ubeta, dsp_data$utau, x_seed)

    # xi testing data
    xi_seed <- 21L
    target_samples_xi <- utest_xi(dsp_data, W, ubeta, phi_specs["mean"], xi_seed)

    # gamma categorical testing data
    gam_cat_seed <- 902L
    out_utest_gamma_categ <- utest_gamma_categ(dsp_data, W, xi, ubeta, gam_cat_seed)

    # collect seeds
    seed_vals <- c(gamma_categ = gam_cat_seed,
                   phi         = phi_seed,
                   W           = w_seed,
                   X           = x_seed,
                   xi          = xi_seed)

    # collect testing objects
    test_data <- list(input_gamma_specs          = out_utest_gamma_categ$input_specs,
                      input_ubeta                = ubeta,
                      input_w                    = W,
                      input_x                    = out_utest_x$test_data,
                      input_xi                   = xi,
                      target_data_gamma_categ    = out_utest_gamma_categ$target_data,
                      target_data_phi            = out_utest_phi$target_data,
                      target_samples_gamma_categ = out_utest_gamma_categ$target_samples,
                      target_samples_phi         = out_utest_phi$target_samples,
                      target_samples_w           = target_samples_w,
                      target_samples_x           = out_utest_x$target_output,
                      target_samples_xi          = target_samples_xi,
                      seed_vals                  = seed_vals,
                      epsilon                    = 1e-12)

    # pass testing data to C++ testing driver
    utest_cpp_(U                 = dsp_data$U,
               X_rcpp            = dsp_data$intercourse$X,
               w_day_blocks      = dsp_data$w_day_blocks,
               w_to_days_idx     = w_to_days_idx_old,
               w_cyc_to_subj_idx = dsp_data$w_cyc_to_subj_idx,
               subj_day_blocks   = dsp_data$subj_day_blocks,
               day_to_subj_idx   = dsp_data$day_to_subj_idx,
               gamma_specs       = gamma_hyper_list,
               phi_specs         = phi_specs,
               x_miss_cyc        = dsp_data$intercourse$miss_cyc,
               x_miss_day        = dsp_data$intercourse$miss_day,
               utau_rcpp         = dsp_data$utau,
               tau_coefs         = dsp_data$tau_fit,
               fw_len            = 5,
               n_burn            = 0,
               n_samp            = n_samp,
               test_data         = test_data)
}
