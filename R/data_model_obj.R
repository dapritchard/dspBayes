get_model_obj <- function(comb_dat, var_nm, dsp_model) {

    # construct lists with each element a vector indexing all the observations
    # in the daily data for a given subject or cycle, respectively
    subj_idx_list <- get_subj_idx(comb_dat, var_nm)
    cyc_idx_list <- get_cyc_idx(comb_dat, var_nm)

    # extract pregnancy information,
    Y <- day_to_cyc_preg(comb_dat, var_nm, cyc_idx_list)

    # extract intercourse data and index the missing values, and then initialize
    # the missing to "no intercourse" status
    X <- comb_dat[[var_nm$sex]]
    miss_x_idx <- is.na(X) %>% which
    X[miss_x_idx] <- FALSE

    U  = expand_model_rhs(comb_dat, dsp_model)
    #### TODO check if data is collinear or constant within outcome ####

    cov_miss_info <- get_missing_var_info(U, dsp_model)

    var_categ_status <- get_var_categ_status(cov_miss_info)

    list(subj_idx_list    = subj_idx_list,
         cyc_idx_list     = cyc_idx_list,
         cov_miss_info    = cov_miss_info,
         miss_x_idx       = miss_x_idx,
         var_categ_status = var_categ_status,
         Y                = Y,
         X                = X,
         U                = U)
}




get_subj_idx <- function(comb_dat, var_nm) {
    id_vec <- unique(comb_dat[[var_nm$id]])
    unique_id_vec <- unique(id_vec)
    lapply(unique_id_vec, function(x) which(id_vec == x))
}




get_cyc_idx <- function(comb_dat, var_nm) {

    keypairs <- get_keypairs(comb_dat, NULL, var_nm)
    cyc_idx_list <- vector("list", NROW(keypairs))

    for (i in seq_along(cyc_idx_list)) {

        curr_id <- keypairs[i, var_nm$id]
        curr_cyc <- keypairs[i, var_nm$cyc]
        cyc_idx_list[[i]] <- (comb_dat[[var_nm$id]] == curr_id &
                              comb_dat[[var_nm$cyc]] == curr_cyc) %>% which
    }

    cyc_idx_list
}




day_to_cyc_preg <- function(comb_dat, var_nm, cyc_idx_list) {

    day_preg_vec <- comb_dat[[var_nm$preg]]
    cyc_preg_vec <- vector("logical", length(cyc_idx_list))

    for (i in seq_along(cyc_idx_list)) {

        curr_idx <- cyc_idx_list[[i]]

        curr_preg_vec <- map_to_bool(day_preg_vec[curr_idx])
        if (length(table(curr_preg_vec)) > 1L) {
            msg <- paste0("inconsistent pregnancy data for subject ", curr_id,
                          "and cycle ", curr_cyc)
            stop(msg, call. = FALSE)
        }

        cyc_preg_vec[i] <- curr_preg_vec[1L]
    }

    # ensure that pregnancy is not constant over cycles
    if (all(preg_vec)) {
        stop("all of the cycles resulted in a pregnancy", call. = FALSE)
    } else if (all(! preg_vec)) {
        stop("none of the cycles resulted in a pregnancy", call. = FALSE)
    }

    preg_vec
}




# expand the model RHS in the presence of missing.  Without missing the call
# could simply be `model.matrix(dsp_model, comb_dat)`.  See
# stackoverflow.com/q/5616210 for details.

expand_model_rhs <- function(comb_dat, dsp_model) {
    model.matrix(dsp_model, model.frame(dsp_model, comb_dat, na.action = na.pass))
}




var_categ_status <- get_var_categ_status(cov_miss_info, n_vars) {

    # the k-th element of `var_categ_status` tracks whether the k-th column in
    # the design matrix is a binary or continuous variable
    var_categ_status <- rep(FALSE, n_vars)

    # each iteration looks up the information for one of the unexpanded
    # variables in the data, and if it is categorical changes the status of the
    # elements in `var_categ_status` corresponding to the expanded desing
    # matrix.
    for (curr_var in cov_miss_info) {

        # if current unexpanded variable is categorical, then get the start and
        # end indices corresponding to the design matrix, and set the
        # corresponding indices in `var_categ_status` to TRUE.  The `+ 1L` is
        # because the indices are stored useing 0-based indexing.
        if (curr_var$categ) {
            curr_idx <- seq(curr_var$col_start, curr_var$col_end, 1L) + 1L
            var_categ_status[curr_idx] <- TRUE
        }
    }

    var_categ_status
}
