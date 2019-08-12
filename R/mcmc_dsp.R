dsp <- function(dsp_data,
                n_samp      = 10000L,
                nBurn       = 5000L,
                nThin       = 1L,
                hypGam      = NULL,
                tuningGam   = NULL,
                hypPhi      = NULL,
                tuningPhi   = 0.3,
                trackProg   = "percent",
                progQuants  = seq(0.1, 1.0, 0.1),
                fw_ar_model = FALSE,
                fw_decay    = NULL,
                mday_priors = NULL,
                noninf_model = FALSE) {

    stopifnot(! (fw_ar_model && is.null(fw_decay)))
    # create dummy values to ensure that we don't pass a null in to the C++
    # function
    if (!fw_ar_model) {
        fw_decay = rep(1, (2 * dsp_data$fw_len) - 1)
    }

    # provide default values for the peak intensity days if necessary, and
    # ensure that the values sum to 1
    if (is.null(mday_priors)) {
        mday_priors <- rep(1 / dsp_data$fw_len, dsp_data$fw_len)
    }
    else if (sum(mday_priors) != 1) {
        mday_priors <- mday_priors / sum(mday_priors)
    }

    # stub functions for gamma and phi specs
    gamma_hyper_list <- get_gamma_specs(dsp_data, fw_ar_model, noninf_model)
    phi_specs <- get_phi_specs()

    # TODO: need to insert a way to add priors for UGen

    # start timer
    start_time <- proc.time()

    out <- dsp_(u_rcpp            = dsp_data$U,
                x_rcpp            = dsp_data$x_vals,
                w_day_blocks      = dsp_data$w_day_blocks,
                w_to_days_idx     = dsp_data$w_to_days_idx,
                w_cyc_to_subj_idx = dsp_data$w_cyc_to_subj_idx,
                subj_day_blocks   = dsp_data$subj_day_blocks,
                day_to_subj_idx   = dsp_data$day_to_subj_idx,
                gamma_specs       = gamma_hyper_list,
                phi_specs         = phi_specs,
                x_miss            = dsp_data$x_miss,
                sex_miss_to_w     = dsp_data$sex_miss_to_w,
                sex_miss_info     = dsp_data$sex_miss_info,
                utau_rcpp         = dsp_data$utau,
                tau_coefs         = dsp_data$tau_fit,
                u_miss_info       = dsp_data$u_miss_info,
                u_miss_type       = dsp_data$u_miss_type,
                u_preg_map        = dsp_data$cov_miss_w_idx,
                u_sex_map         = dsp_data$cov_miss_x_idx,
                fw_len            = dsp_data$fw_len,
                n_burn            = 0L,
                n_samp            = n_samp,
                fw_decay          = fw_decay,
                log_mday_priors   = log(mday_priors))

    # end timer
    run_time <- proc.time() - start_time

    # transpose data
    coefs_trans <- matrix(out$coefs,
                          nrow      = n_samp,
                          ncol      = ncol(dsp_data$U),
                          byrow     = TRUE,
                          dimnames  = list(NULL, colnames(dsp_data$U)))
    xi_trans <- matrix(out$xi,
                       nrow   = n_samp,
                       byrow  = TRUE)

    list(coefs    = coefs_trans,
         xi       = xi_trans,
         phi      = out$phi,
         ugen     = out$ugen,
         mday     = out$mday,
         mu       = out$mu,
         nu       = out$nu,
         delta    = out$delta,
         run_time = run_time)
}




# FIXME: remove this!
dsp_v2 <- function(dsp_data,
                n_samp      = 10000L,
                nBurn       = 5000L,
                nThin       = 1L,
                hypGam      = NULL,
                tuningGam   = NULL,
                hypPhi      = NULL,
                tuningPhi   = 0.3,
                trackProg   = "percent",
                progQuants  = seq(0.1, 1.0, 0.1),
                fw_ar_model = FALSE) {

    # stub functions for gamma and phi specs
    gamma_hyper_list <- get_gamma_specs(dsp_data, fw_ar_model)
    for (i in 1:5) gamma_hyper_list[[i]]["type"] <- 0
    phi_specs <- get_phi_specs()

    # TODO: need to insert a way to add priors for UGen

    # start timer
    start_time <- proc.time()

    out <- dsp_(u_rcpp            = dsp_data$U,
                x_rcpp            = dsp_data$x_vals,
                w_day_blocks      = dsp_data$w_day_blocks,
                w_to_days_idx     = dsp_data$w_to_days_idx,
                w_cyc_to_subj_idx = dsp_data$w_cyc_to_subj_idx,
                subj_day_blocks   = dsp_data$subj_day_blocks,
                day_to_subj_idx   = dsp_data$day_to_subj_idx,
                gamma_specs       = gamma_hyper_list,
                phi_specs         = phi_specs,
                x_miss            = dsp_data$x_miss,
                sex_miss_to_w     = dsp_data$sex_miss_to_w,
                sex_miss_info     = dsp_data$sex_miss_info,
                utau_rcpp         = dsp_data$utau,
                tau_coefs         = dsp_data$tau_fit,
                u_miss_info       = dsp_data$u_miss_info,
                u_miss_type       = dsp_data$u_miss_type,
                u_preg_map        = dsp_data$cov_miss_w_idx,
                u_sex_map         = dsp_data$cov_miss_x_idx,
                fw_len            = dsp_data$fw_len,
                n_burn            = 0L,
                n_samp            = n_samp)

    # end timer
    run_time <- proc.time() - start_time

    # transpose data
    coefs_trans <- matrix(out$coefs,
                          nrow      = n_samp,
                          ncol      = ncol(dsp_data$U),
                          byrow     = TRUE,
                          dimnames  = list(NULL, colnames(dsp_data$U)))
    xi_trans <- matrix(out$xi,
                       nrow   = n_samp,
                       byrow  = TRUE)

    list(coefs    = coefs_trans,
         xi       = xi_trans,
         phi      = out$phi,
         ugen     = out$ugen,
         mu       = out$mu,
         nu       = out$nu,
         delta    = out$delta,
         run_time = run_time)
}
