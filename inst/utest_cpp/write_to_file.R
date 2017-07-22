write_to_file <- function(dsp_data, dir_nm) {

    # construct directory into which we will place data
    path_above_nm <- file.path("inst", "utest_cpp")
    path_nm <- file.path("inst", "utest_cpp", dir_nm)

    # create data directory.  Throw exception if (i) we can't find the level
    # above the directory or (ii) if the directory already exists (so that we
    # don't risk clobbering existing data)
    if (! file.exists(path_above_nm)) {
        stop(paste0(path_above_nm, " relative directory doesn't exist"))
    } else if (file.exists(path_nm)) {
        stop(paste0("file ", dir_nm, " already exists, aborting"))
    } else {
        system(paste0("mkdir ", path_nm))
    }


    fw_len <- 5
    n_burn <- 0
    n_samp <- 10

    writeBin(as.vector(dsp_data$U), file.path(path_nm, "input_U"))

    writeBin(dsp_data$X, file.path(path_nm, "input_X_rcpp"))

    # preg_cyc_beg_idx <- sapply(dsp_data$w_day_blocks, function(x) x["beg_idx"])
    # preg_cyc_n_days <- sapply(dsp_data$w_day_blocks, function(x) x["n_days"])
    # preg_cyc_subj_idx <- sapply(dsp_data$w_day_blocks, function(x) x["subj_idx"])
    # writeBin(preg_cyc_beg_idx, file.path(path_nm, "input_w_day_blocks_beg"))
    # writeBin(preg_cyc_n_days, file.path(path_nm, "input_w_day_blocks_n"))
    # writeBin(preg_cyc_subj_idx, file.path(path_nm, "input_w_day_blocks_subj"))
    writeBin(unlist(dsp_data$w_day_blocks), file.path(path_nm, "input_w_day_blocks"))

    writeBin(dsp_data$w_to_days_idx, file.path(path_nm, "input_w_to_days_idx"))

    writeBin(dsp_data$w_cyc_to_cyc_idx, file.path(path_nm, "input_w_cyc_to_cyc_idx"))

    # subj_beg_idx <- sapply(dsp_data$subj_day_blocks, function(x) x["beg_idx"])
    # subj_n_days <- sapply(dsp_data$subj_day_blocks, function(x) x["n_days"])
    # writeBin(subj_beg_idx, file.path(path_nm, "input_subj_day_blocks_beg"))
    # writeBin(subj_n_days, file.path(path_nm, "input_subj_day_blocks_n"))
    writeBin(unlist(dsp_data$subj_day_blocks), file.path(path_nm, "input_subj_day_blocks"))

    writeBin(dsp_data$day_to_subj_idx, file.path(path_nm, "input_day_to_subj_idx"))

    gamma_specs <- get_gamma_specs(dsp_data) %>% unlist
    writeBin(gamma_specs, file.path(path_nm, "input_gamma_specs"))

    phi_specs <- get_phi_specs(dsp_data)
    writeBin(phi_specs, file.path(path_nm, "input_phi_specs"))

    writeBin(c(fw_len, n_burn, n_samp), file.path(path_nm, "input_misc_specs"))
}
