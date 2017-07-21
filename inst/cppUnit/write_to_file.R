write_to_file <- function(dsp_data, dir_nm) {

    # construct directory into which we will place data
    path_nm <- file.path("inst", "cppUnit", dir_nm)

    # create data directory.  Throw exception if directory exists so that we
    # don't risk clobbering existing data, or if we can't find the directory
    if (! file.exists(file.path("inst", "cppUnit"))) {
        stop("inst/cppUnit relative directory doesn't exist")
    } else if (file.exists(path_nm)) {
        stop(paste0("file ", dir_nm, " already exists, aborting"))
    } else {
        system(paste0("mkdir inst/cppUnit/", dir_nm))
    }


    writeBin(as.vector(dsp_data$U), file.path(path_nm, "input_U.dat"))

    writeBin(X, file.path(path_nm, "input_X.dat"))

    preg_cyc_beg_idx <- sapply(dsp_data$preg_cyc_list, function(x) x["beg_idx"])
    preg_cyc_n_days <- sapply(dsp_data$preg_cyc_list, function(x) x["n_days"])
    preg_cyc_subj_idx <- sapply(dsp_data$preg_cyc_list, function(x) x["subj_idx"])
    writeBin(preg_cyc_beg_idx, file.path(path_nm, "input_preg_cyc_beg.dat"))
    writeBin(preg_cyc_n_days, file.path(path_nm, "input_preg_cyc_n.dat"))
    writeBin(preg_cyc_subj_idx, file.path(path_nm, "input_preg_cyc_subj.dat"))

    writeBin(dsp_data$w_days_idx, file.path(path_nm, "input_w_days_idx.dat"))

    writeBin(dsp_data$w_cyc_idx, file.path(path_nm, "input_w_cyc_idx.dat"))

    subj_beg_idx <- sapply(dsp_data$subj_idx_list, function(x) x["beg_idx"])
    subj_n_days <- sapply(dsp_data$subj_idx_list, function(x) x["n_days"])
    writeBin(subj_beg_idx, file.path(path_nm, "input_subj_beg.dat"))
    writeBin(subj_n_days, file.path(path_nm, "input_subj_n.dat"))

    # TODO: gamma specs

    # TODO: phi hyper

    # TODO: fw_len, n_burn, n_samp, n_thin
}





# x <- rnorm(5)

# writeBin(x, "inst/cppUnit/test_double.dat")
