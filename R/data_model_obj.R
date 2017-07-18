derive_model_obj <- function(comb_dat, var_nm, dsp_model) {


    preg_cyc_list <- get_cyc_list(comb_dat, var_nm)
    preg_days_idx <- get_w_days_idx(preg_cyc_list)
    preg_cyc_idx <- get_cyc_idx(preg_cyc_list)

    subj_idx_list <- get_subj_idx(comb_dat, var_nm)
    days_to_subj_idx <- get_days_to_subj(subj_idx_list)

    # extract intercourse data and index the missing values, and then initialize
    # the missing to "no intercourse" status (using 0/1 coding for no/yes).
    # Note that the expression for `miss_x_idx` does the right thing when there
    # are no missing, which is to provide an `integer(0)` vector.
    miss_x_bool <- is.na(comb_dat[[var_nm$sex]])
    X <- as.integer(! miss_x_bool)
    miss_x_idx <- which(miss_x_bool) - 1L

    U  = expand_model_rhs(comb_dat, dsp_model)
    #### TODO check if data is collinear or constant within outcome ####

    cov_miss_info <- get_missing_var_info(U, dsp_model)

    var_categ_status <- get_var_categ_status(cov_miss_info)

    list(preg_cyc_list    = preg_cyc_list,
         preg_days_idx    = preg_days_idx,
         preg_cyc_idx     = preg_cyc_idx,
         subj_idx_list    = subj_idx_list,
         days_to_subj_idx = days_to_subj_idx,
         miss_x_idx       = miss_x_idx,
         cov_miss_info    = cov_miss_info,
         var_categ_status = var_categ_status,
         X                = X,
         U                = U)
}




get_subj_idx <- function(comb_dat, var_nm) {
    id_vec <- comb_dat[[var_nm$id]]
    unique_id_vec <- unique(id_vec)
    lapply(unique_id_vec, function(x) {
        curr_idx <- which(id_vec == x)
        # subtract 1 to convert to 0-based indexing
        c(beg_idx = head(curr_idx, 1L) - 1L, n_days = length(curr_idx))
    })
}




get_days_to_subj <- function(subj_idx_list) {
    n_days <- sapply(subj_idx_list, function(x) x["n_days"])
    rep(seq_along(subj_idx_list), n_days)
}




get_cyc_list <- function(comb_dat, var_nm) {

    keypairs <- get_keypairs(comb_dat, NULL, var_nm)
    cyc_idx_list <- vector("list", NROW(keypairs))
    ctr <- 1L

    for (i in seq_along(cyc_idx_list)) {

        curr_id <- keypairs[i, var_nm$id]
        curr_cyc <- keypairs[i, var_nm$cyc]
        curr_idx <- which(comb_dat[[var_nm$id]] == curr_id &
                          comb_dat[[var_nm$cyc]] == curr_cyc)

        # obtain whether pregnancy occured during the current cycle and store as
        # `curr_preg_bool`
        curr_preg_vec <- map_vec_to_bool(comb_dat[curr_idx, var_nm$preg])
        if (length(table(curr_preg_vec)) > 1L) {
            msg <- paste0("inconsistent pregnancy data for subject ", curr_id,
                          "and cycle ", curr_cyc)
            stop(msg, call. = FALSE)
        }

        # case: a pregnancy occurred so record the cycle information into
        # `cyc_idx_list`
        if (curr_preg_vec[1L]) {

            # obtain the
            beg_idx <- head(curr_idx, 1L)
            n_days <- length(curr_idx)
            subj_idx <- which(unique(comb_dat[[var_nm$id]]) == curr_id)

            # subtract 1 to convert to 0-based indexing
            cyc_idx_list[[ctr]] <- c(beg_idx  = beg_idx - 1L,
                                     n_days   = n_days,
                                     subj_idx = subj_idx - 1L,
                                     cyc_idx  = i - 1L)
            ctr <- ctr + 1L
        }

    } # end step through all cycles loop

    cyc_idx_list[1:(ctr - 1L)]
}




get_w_days_idx <- function(preg_cyc_list) {
    idx_list <- lapply(preg_cyc_list, function(x) {
        x["beg_idx"] : (x["beg_idx"] + x["n_days"] - 1L)
    })
    unlist(idx_list)
}




get_cyc_idx <- function(preg_cyc_list) {

    cyc_idx <- vector("integer", length(preg_cyc_list))
    for (i in seq_along(preg_cyc_list)) {
        cyc_idx[i] <- preg_cyc_list[[i]]["cyc_idx"]
    }

    cyc_idx
}




# day_to_cyc_preg <- function(comb_dat, var_nm, cyc_idx_list) {

#     day_preg_vec <- comb_dat[[var_nm$preg]]
#     cyc_preg_vec <- vector("logical", length(cyc_idx_list))

#     for (i in seq_along(cyc_idx_list)) {

#         curr_idx <- cyc_idx_list[[i]]

#         curr_preg_vec <- map_vec_to_bool(day_preg_vec[curr_idx])
#         if (length(table(curr_preg_vec)) > 1L) {
#             msg <- paste0("inconsistent pregnancy data for subject ", curr_id,
#                           "and cycle ", curr_cyc)
#             stop(msg, call. = FALSE)
#         }

#         cyc_preg_vec[i] <- curr_preg_vec[1L]
#     }

#     # # ensure that pregnancy is not constant over cycles
#     # if (all(cyc_preg_vec)) {
#     #     stop("all of the cycles resulted in a pregnancy", call. = FALSE)
#     # } else if (all(! cyc_preg_vec)) {
#     #     stop("none of the cycles resulted in a pregnancy", call. = FALSE)
#     # }

#     cyc_preg_vec
# }




# expand the model RHS in the presence of missing.  Without missing the call
# could simply be `model.matrix(dsp_model, comb_dat)`.  See
# stackoverflow.com/q/5616210 for details.

expand_model_rhs <- function(comb_dat, dsp_model) {
    model.matrix(dsp_model, model.frame(dsp_model, comb_dat, na.action = na.pass))
}




get_var_categ_status <- function(cov_miss_info, n_vars) {

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
