get_gamma_specs <- function(dsp_data, fw_ar_model, noninf_model) {

    stopifnot(!is.null(dsp_data$fw_len))
    mh_delta_factor <- 4
    if (dsp_data$fw_len == 5L) {
        prior_expected_val <- c(0.115, 0.196, 0.288, 0.340, 0.264)
        mh_delta_vals <- c(0.4, 0.2, 0.2, 0.2, 0.4) * mh_delta_factor
    }
    else if (dsp_data$fw_len == 7L) {
        prior_expected_val <- c(0.059, 0.115, 0.196, 0.288, 0.340, 0.264, 0.104)
        mh_delta_vals <- c(0.4, 0.4, 0.2, 0.2, 0.2, 0.4, 0.4) * mh_delta_factor
    }
    else if (dsp_data$fw_len == 9L) {
        prior_expected_val <- c(0.026, 0.059, 0.115, 0.196, 0.288, 0.340, 0.264, 0.104, 0.017)
        mh_delta_vals <- c(0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4) * mh_delta_factor
    }
    else {
        stop("haven't defined ", dsp_data$fw_len, " number of FW days")
    }
    val_factor <- 10

    # mean is `prior_expected_val`, and variance is `prior_expected_val / val_factor`
    if (noninf_model) {
        prior_expected_val <- rep(0.25, dsp_data$fw_len)
        val_factor <- 5
    }
    a_lookup <- prior_expected_val * val_factor
    b_lookup <- rep(val_factor, length(prior_expected_val))

    # some temporary glue code.  create gamma specs
    gamma_hyper_list <- vector("list", ncol(dsp_data$U))
    for (i in seq_len(ncol(dsp_data$U))) {
        gamma_hyper_list[[i]] <- c(type     = 1,
                                   h        = i - 1,
                                   hyp_a    = a_lookup[i],
                                   hyp_b    = b_lookup[i],
                                   hyp_p    = 0,
                                   # hyp_a    = 1,
                                   # hyp_b    = 1,
                                   # hyp_p    = 0.5,
                                   bnd_l    = 0,
                                   bnd_u    = Inf,
                                   mh_p     = 0.1,
                                   mh_delta = c(mh_delta_vals, (rep(0.5, 20)))[i])
    }

    if (fw_ar_model) {       # FIXME: hardcoded number of FW days!!
        for (i in seq_along(prior_expected_val)) {
            gamma_hyper_list[[i]]["type"] <- 3
            # hyp_p <- 0
        }
    }

    gamma_hyper_list
}


get_phi_specs <- function(dsp_data) {
    c(c1 = 1, c2 = 1, delta = 0.1, mean = 1)
}
