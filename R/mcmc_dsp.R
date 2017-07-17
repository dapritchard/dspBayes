dsp <- function(dspDat,
                nSamp      = 1e4,
                nBurn      = 5e3,
                nThin      = 1L,
                hypGam     = NULL,
                tuningGam  = NULL,
                hypPhi     = NULL,
                tuningPhi  = 0.3,
                trackProg  = "percent",
                progQuants = seq(0.1, 1.0, 0.1)) {


    dsp_sampler(U,
                preg_cyc_list,
                w_days_idx,
                w_cyc_idx,
                fw_len,
                xi_initial,
                subj_days_list,
                gamma_specs_list,
                phi_hyper_list)
}
