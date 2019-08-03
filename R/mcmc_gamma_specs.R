get_gamma_specs <- function(dsp_data, fw_ar_model) {

    true_vals <- c(0.115, 0.196, 0.288, 0.34, 0.264)
    a_lookup <- true_vals^2
    b_lookup <- true_vals

    # some temporary glue code.  create gamma specs
    gamma_hyper_list <- vector("list", ncol(dsp_data$U))
    for (i in seq_len(ncol(dsp_data$U))) {
        stopifnot(i <= 5L)
        gamma_hyper_list[[i]] <- c(type     = 1,
                                   h        = i - 1,
                                   hyp_a    = a_lookup[i],
                                   hyp_b    = b_lookup[i],
                                   hyp_p    = 0.00,
                                   # hyp_a    = 1,
                                   # hyp_b    = 1,
                                   # hyp_p    = 0.5,
                                   bnd_l    = 0,
                                   bnd_u    = Inf,
                                   mh_p     = 0.1,
                                   mh_delta = c(c(0.4, 0.2, 0.2, 0.2, 0.4)/0.25, (rep(0.5, 20)))[i])
    }

    if (fw_ar_model) {       # FIXME: hardcoded number of FW days!!
        for (i in 1L:5L) {
            gamma_hyper_list[[i]]["type"] <- 3
            hyp_p <- 0
        }
    }

    gamma_hyper_list
}


get_phi_specs <- function(dsp_data) {
    c(c1 = 1, c2 = 1, delta = 0.1, mean = 1)
}
