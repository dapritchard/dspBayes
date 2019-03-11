get_gamma_specs <- function(dsp_data, fw_ar_model) {

    # some temporary glue code.  create gamma specs
    gamma_hyper_list <- vector("list", ncol(dsp_data$U))
    for (i in seq_len(ncol(dsp_data$U))) {
        gamma_hyper_list[[i]] <- c(type     = 0,
                                   h        = i - 1,
                                   hyp_a    = 1,
                                   hyp_b    = 1,
                                   hyp_p    = 0.5,
                                   bnd_l    = 0,
                                   bnd_u    = Inf,
                                   mh_p     = 0.1,
                                   mh_delta = 0.1)
    }

    if (fw_ar_model) {       # FIXME: hardcoded number of FW days!!
        for (i in 1L:5L) {
            gamma_hyper_list[[i]]["type"] <- 3
        }
    }

    gamma_hyper_list
}


get_phi_specs <- function(dsp_data) {
    c(c1 = 1, c2 = 1, delta = 0.1, mean = 1)
}
