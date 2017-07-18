dsp <- function(dsp_data,
                nSamp      = 1e4,
                nBurn      = 5e3,
                nThin      = 1L,
                hypGam     = NULL,
                tuningGam  = NULL,
                hypPhi     = NULL,
                tuningPhi  = 0.3,
                trackProg  = "percent",
                progQuants = seq(0.1, 1.0, 0.1)) {

    # some scratch glue code.  create gamma specs
    gamma_hyper_list <- vector("list", ncol(dsp_data$U))
    for (i in seq_len(ncol(dsp_data$U))) {
        gamma_hyper_list[[i]] <- c(type = 0,
                                   h = i,
                                   hyp_a = 1,
                                   hyp_b = 1,
                                   hyp_p = 0.5,
                                   bnd_l = 0,
                                   bnd_u = Inf)
    }
    # some scratch glue code.  create phi specs
    phi_hyper <- c(c1 = 1, c2 = 1, delta = 0.1, mean = 1)



    dsp_sampler(U           = dsp_data$U,
                X           = dsp_data$X,
                preg_cyc    = dsp_data$preg_cyc_list,
                w_days_idx  = dsp_data$preg_days_idx,
                w_cyc_idx   = dsp_data$preg_cyc_idx,
                subj_days   = dsp_data$subj_idx_list,
                subj_idx    = dsp_data$days_to_subj_idx,
                gamma_specs = gamma_hyper_list,
                phi_hyper   = phi_hyper,
                fw_len      = 5,
                n_burn      = 0,
                n_samp      = 5)

}
