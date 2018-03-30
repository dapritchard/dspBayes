derive_model_obj <- function(filtered_dat_list,
                             var_nm,
                             fw_incl,
                             dsp_model,
                             use_na,
                             tau_fit) {

    # bind `comb_dat` to the combined data stored in
    comb_dat <- filtered_dat_list$comb_dat

    # miss_sex_cyc_idx <- filtered_dat$miss_sex_cyc_idx
    # day_zero_sex     <- filtered_dat$day_zero_sex

    w_day_blocks <- get_w_day_blocks(comb_dat, var_nm)
    w_to_days_idx <- get_w_to_days_idx(comb_dat, var_nm)
    w_cyc_to_subj_idx <- get_w_cyc_to_subj_idx(w_day_blocks)

    subj_day_blocks <- get_subj_day_blocks(comb_dat, var_nm)
    day_to_subj_idx <- get_day_to_subj_idx(subj_day_blocks)

    U <- expand_model_rhs(comb_dat, dsp_model)
    #### TODO check if data is collinear or constant within outcome ####

    # obtain intercourse values with missing filled in, and track which values
    # were originally missing
    x_vals <- get_intercourse_data_v2(comb_dat, var_nm, fw_incl)
    x_miss <- is.na(comb_dat[[var_nm$sex]]) %>% as.integer

    # a vector with length equal to the number of missing intercourse days,
    # where each day is the index in W if the day corresponds to a pregnancy
    # cycle, or -1 otherwise
    sex_miss_to_w <- get_sex_miss_to_w(comb_dat, var_nm)


    sex_miss_info <- get_sex_info(comb_dat, var_nm,  use_na)

    # # set the indices in `xmiss` now that non-sex days have been removed
    # xmiss <- map_xmiss_to_x(comb_dat, var_nm, xmiss)

    cov_col_miss_info <- get_cov_col_miss_info(U, dsp_model, use_na)
    cov_row_miss_info <- get_cov_row_miss_info(comb_dat, var_nm, U, cov_col_miss_info)
    u_miss_info <- get_u_miss_info(cov_col_miss_info, cov_row_miss_info)
    u_miss_type <- get_var_categ_status(cov_col_miss_info)
    u_miss_filled_in <- get_u_miss_filled_in(U, cov_col_miss_info)

    # note that we have to perform this after `U` has missing values filled in
    utau <- get_utau(u_miss_filled_in, tau_fit, use_na)

    list(w_day_blocks      = w_day_blocks,
         w_to_days_idx     = w_to_days_idx,
         w_cyc_to_subj_idx = w_cyc_to_subj_idx,
         subj_day_blocks   = subj_day_blocks,
         day_to_subj_idx   = day_to_subj_idx,
         x_vals            = x_vals,
         x_miss            = x_miss,
         sex_miss_to_w     = sex_miss_to_w,
         cov_col_miss_info = cov_col_miss_info,
         tau_fit           = tau_fit,
         utau              = utau,
         u_miss_info       = u_miss_info,
         u_miss_type       = u_miss_type,
         cov_miss_w_idx    = cov_row_miss_info$cov_miss_w_idx,
         cov_miss_x_idx    = cov_row_miss_info$cov_miss_x_idx,
         U                 = u_miss_filled_in)
}




get_w_day_blocks <- function(comb_dat, var_nm) {

    keypairs <- get_keypairs(comb_dat, NULL, var_nm)
    cyc_idx_list <- vector("list", NROW(keypairs))

    # each iteration checks a block of days in `comb_dat` corresponding to the
    # current cycles to see if a pregnancy occured during the cycle.  If one did
    # occur, then the info for the block is recorded in the `ctr`-th element of
    # `cyc_idx_list`.
    ctr <- 1L
    for (i in seq_along(cyc_idx_list)) {

        # TODO: modify this using `get_cycle_idx`

        # get the indices in `comb_dat` corresponding to the current cycle
        curr_id <- keypairs[i, var_nm$id]
        curr_cyc <- keypairs[i, var_nm$cyc]
        curr_idx <- which(comb_dat[[var_nm$id]] == curr_id &
                          comb_dat[[var_nm$cyc]] == curr_cyc)

        # obtain whether pregnancy occured during the current cycle and store as
        # `curr_preg_bool`
        curr_preg_vec <- map_vec_to_bool(comb_dat[curr_idx, var_nm$preg])
        if (length(table(curr_preg_vec)) > 1L) {
            stop("inconsistent pregnancy data for subject ", curr_id,
                 "and cycle ", curr_cyc, call. = FALSE)
        }

        # case: a pregnancy occurred so record the cycle information into
        # `cyc_idx_list`
        if (curr_preg_vec[1L]) {

            # obtain the needed values for the current block
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




# get_w_to_days_idx <- function(preg_cyc_list) {
#     idx_list <- lapply(preg_cyc_list, function(x) {
#         x["beg_idx"] : (x["beg_idx"] + x["n_days"])
#     })
#     # now back to 0-based indexing
#     unlist(idx_list)
# }




# calculates the (0-based) day-specific indices for the days in which a
# pregnancy occured during the corresponding cycle and adds a sentinal value at
# the end of the vector used to signal when there are no more days occuring in a
# pregnancy cycle
#
# PRE: `comb_dat` is a data frame with one of the columns having the name given
# by `var_nm$preg`.  It is assumed that there is no missing in
# `comb_dat[[var_nm$preg]]`, although these would be removed anyway.

get_w_to_days_idx <- function(comb_dat, var_nm) {
    comb_dat[[var_nm$preg]] %>% map_vec_to_bool %>% which %>% c(., 0L) %>% `-`(., 1L)
}




get_w_cyc_to_subj_idx <- function(preg_cyc_list) {
    sapply(preg_cyc_list, function(x) x["subj_idx"]) %>% structure(., names = NULL)
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




get_subj_day_blocks <- function(comb_dat, var_nm) {

    lapply(unique( comb_dat[[var_nm$id]] ), function(x) {

        curr_idx <- which(comb_dat[[var_nm$id]] == x)

        # subtract 1 to convert to 0-based indexing
        c(beg_idx = head(curr_idx, 1L) - 1L,
          n_days  = length(curr_idx))
    })
}




get_day_to_subj_idx <- function(subj_day_blocks) {
    n_days <- sapply(subj_day_blocks, function(x) x["n_days"])
    # subtract 1 to convert to 0-based indexing
    rep(seq_along(subj_day_blocks), n_days) - 1L
}




# expand the model RHS in the presence of missing.  Without missing the call
# could simply be `model.matrix(dsp_model, comb_dat)`.  See
# stackoverflow.com/q/5616210 for details.

expand_model_rhs <- function(comb_dat, dsp_model) {
    model.matrix(dsp_model, model.frame(dsp_model, comb_dat, na.action = na.pass))
}



# # TODO: what is this function? should it be used?
# get_var_categ_status <- function(cov_miss_info, n_vars) {

#     # the k-th element of `var_categ_status` tracks whether the k-th column in
#     # the design matrix is a binary or continuous variable
#     var_categ_status <- rep(FALSE, n_vars)

#     # each iteration looks up the information for one of the unexpanded
#     # variables in the data, and if it is categorical changes the status of the
#     # elements in `var_categ_status` corresponding to the expanded design
#     # matrix.
#     for (curr_var in cov_miss_info) {

#         # if current unexpanded variable is categorical, then get the start and
#         # end indices corresponding to the design matrix, and set the
#         # corresponding indices in `var_categ_status` to TRUE.  The `+ 1L` is
#         # because the indices are stored using 0-based indexing.
#         if (curr_var$categ) {
#             curr_idx <- seq(curr_var$col_start, curr_var$col_end, 1L) + 1L
#             var_categ_status[curr_idx] <- TRUE
#         }
#     }

#     var_categ_status
# }


get_var_categ_status <- function(cov_col_miss_info) {
    if (length(cov_col_miss_info) > 0L) {
        sapply(cov_col_miss_info, function(x) x["categ"] %>% as.integer)
    } else {
        vector("integer", 0L)
    }
}




# TODO: no longer need the `miss_cyc` and `miss_day` portions of this code.  The
# only thing we need is the portion that initializes missing intercourse values.

# TODO: no longer use this at all?

get_intercourse_data <- function(comb_dat, var_nm, fw_incl) {

    # a list with each element a vector of the days-specific indices
    # corresponding to one of the cycles
    cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)

    # store intercourse data as a binary variable.  Missingness is preserved.
    X <- map_vec_to_bool(comb_dat[, var_nm$sex]) %>% as.integer
    # TODO: use use_na
    # TODO: check if there are any cycles with a pregnancy and only 1 missing
    # day, and all other days are non-intercourse.  In this case X_ijk must be
    # an intercourse.
    x_miss_bool <- is.na(X)
    x_miss_idx <- which(x_miss_bool)

    # convert sex yesterday to a binary variable and map missings value to the
    # corresponding flags.  Returns NULL if no column "sex_yester" exists
    # (i.e. we are not imputing missing intercourse data)
    sex_yester <- get_sex_yester_coding(comb_dat, var_nm)

    # map the days to pregnancy days.  If the value is 0 then this signals that
    # the day was not a day that occured during a pregnancy cycle.
    preg_day_bool <- map_vec_to_bool( comb_dat[[var_nm$preg]] )
    preg_day_map <- vector("integer", length(X))
    preg_day_map[preg_day_bool] <- seq_len(sum(preg_day_bool))

    # containers to store missing intercourse information
    x_miss_cyc <- vector("list", length(cyc_idx_list))
    x_miss_day <- vector("list", sum(x_miss_bool))

    # maps a missing intercourse day to the index of `x_miss_idx`
    idx_to_x_miss_idx <- vector("integer", length(X))
    idx_to_x_miss_idx[x_miss_bool] <- seq_along(x_miss_idx)

    # map subjects to indices
    id_map <- get_id_map(comb_dat[[var_nm$id]])

    # each iteration checks the current cycle for missing.  If some exists, then
    # add an entry to `x_miss_cyc` detailing the missingness, and fill in the
    # missing elements of X.  Also the maximum amount of missing in a cycle is
    # tracked and stored in `x_n_max_miss`.
    ctr <- 1L
    for (curr_cyc_idx in cyc_idx_list) {

        # which among the current cycle indices are missing, and how many
        curr_miss_bool <- x_miss_bool[curr_cyc_idx]
        curr_miss_idx <- curr_cyc_idx[curr_miss_bool]
        map_to_x_miss_idx <- idx_to_x_miss_idx[curr_miss_idx]
        curr_n_miss <- length(curr_miss_idx)

        # case: at least one day in the current cycle has missing intercourse
        # data
        if (curr_n_miss > 0L) {

            # provide cycle missing information in `x_miss_cyc`.  Subtract 1 to
            # adjust for 0-based indexing.
            x_miss_cyc[[ctr]] <- c(beg_idx  = map_to_x_miss_idx[1L] - 1L,
                                   n_days   = curr_n_miss,
                                   subj_idx = id_map[curr_cyc_idx[1L]] - 1L,
                                   preg_idx = preg_day_map[curr_cyc_idx[1L]] - 1L)

            # set the first missing X in the cycle to 1, and the remaining
            # missing to 0.  The logic for this is that if a pregnancy has
            # occurred in the cycle, then there must be at least 1 day with
            # intercourse, and this guarantees that.
            X[curr_miss_idx[1L]] <- 1L
            X[curr_miss_idx[-1L]] <- 0L

            ctr <- ctr + 1L
        }
    }

    # remove entries corresponding to cycles without any missing
    x_miss_cyc <- x_miss_cyc[seq_len(ctr - 1L)]

    # obtain indices and previous day intercourse status for each missing
    # intercourse observation
    for (i in seq_along(x_miss_idx)) {
        curr_idx <- x_miss_idx[i]
        x_miss_day[[i]] <- c(idx  = curr_idx - 1L,
                             prev = sex_yester[curr_idx])
    }

    # if (length(intercourse_data$miss_day > 0L)) {
    #     tau_data <- get_tau_data()
    # } else {
    #     tau_data <- 1
    # }

    # return intercourse information
    list(X          = X,
         miss_cyc   = x_miss_cyc,
         miss_day   = x_miss_day)
}




get_intercourse_data_v2 <- function(comb_dat, var_nm, fw_incl) {

    # a list with each element a vector of the days-specific indices
    # corresponding to one of the cycles
    cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)

    # store intercourse data as a binary variable.  Missingness is preserved.
    x_vals <- map_vec_to_bool(comb_dat[, var_nm$sex]) %>% as.integer
    # TODO: use use_na
    # TODO: check if there are any cycles with a pregnancy and only 1 missing
    # day, and all other days are non-intercourse.  In this case X_ijk must be
    # an intercourse.

    x_vals[is.na(x_vals)] <- 1L
    x_vals
}




# returns an integer vector with elements that provide indices for the W data in
# pregnancy cycles, or -1L otherwise.  In other words, if the k-th element of
# the return vector has a value of t >= 0, then this means that the k-th missing
# intercourse day corresponds to the t-th day in the W data.  A value of -1L
# means that the k-th missing day corresponds to a non-pregnancy cycle.

get_sex_miss_to_w <- function(comb_dat, var_nm) {

    # map the days to pregnancy days.  If the value is 0 then this signals that
    # the day was not a day that occured during a pregnancy cycle.
    preg_day_bool <- map_vec_to_bool( comb_dat[[var_nm$preg]] )
    preg_day_map <- vector("integer", nrow(comb_dat))
    preg_day_map[preg_day_bool] <- seq_len(sum(preg_day_bool))

    # map missing intercourse days to W days indices.  Subtract 1 for 0-based
    # indexing.  Days that took place in non-pregnancy cycles have a value of
    # -1.
    sex_miss_bool <- is.na(comb_dat[[var_nm$sex]])
    sex_miss_to_w <- preg_day_map[sex_miss_bool] - 1L

    sex_miss_to_w
}




get_subj_idx_list <- function(dataset, var_nm) {

    # container to store the indices in.  Set the length to the maximum possible
    # elements that it could need.
    n <- NROW(dataset)
    subj_idx_list <- vector("list", n)

    # bind variable name to the data for convenience
    id_vec <- dataset[[var_nm$id]]

    # tracks the number of cycles in the data, and the index of the first
    # observation in the current cycle
    ctr <- 1L
    start <- 1L

    # each iteration steps through the current id/cycle keypair until a new
    # cycle is found, and then saves a vector of the indices for the cycle to
    # `subj_idx_list`
    while (start <= n) {

        # obtain the id for the current cycle
        curr_id <- id_vec[start]

        # `end` is used to search for one past the last day for the subject
        end <- start + 1L

        # increment `end` until we get one past the last day for the subject
        while ((end <= n) && (id_vec[end] == curr_id)) {
            end <- end + 1L
        }

        # the current subject index are the indices between `start` and one
        # before `end`.  Updatate `start` and `ctr`.
        subj_idx_list[[ctr]] <- start : (end - 1L)
        start <- end
        ctr <- ctr + 1L
    }

    # return the data shortened to the number of cycles
    subj_idx_list[1:(ctr - 1L)]
}




# returns a list such that each element is a vector providing the indices in
# `dataset` for a given cycle.  The cycles are in order within a vector and
# across vectors, so that performing `unlist` on the return object yields a
# vector with values `1, 2, ..., nrow(dataset)`.
#
# PRE: assumes `dataset` is a data.frame with columns for id and cycle and with
# names as given in `var_nm`.  The function strongly relies on the fact that
# `dataset` is already sorted on the id/cycle keypairs.

# TODO: remove?

get_cyc_idx_list_old <- function(dataset, var_nm) {

    # container to store the indices in.  Set the length to the maximum possible
    # elements that it could need.
    n <- NROW(dataset)
    cyc_idx_list <- vector("list", n)

    # bind variable names to the data for convenience
    id_vec <- dataset[[var_nm$id]]
    cyc_vec <- dataset[[var_nm$cyc]]

    # tracks the number of cycles in the data, and the index of the first
    # observation in the current cycle
    ctr <- 1L
    start <- 1L

    # each iteration steps through the current id/cycle keypair until a new
    # cycle is found, and then saves a vector of the indices for the cycle to
    # `cyc_idx_list`
    while (start <= n) {

        # obtain the id/cycle keypair for the current cycle
        curr_id <- id_vec[start]
        curr_cyc <- cyc_vec[start]

        # `end` is used to search for one past the last day in the cycle
        end <- start + 1L

        # increment `end` until we get one past the last day in the cycle
        while ((end <= n) && (cyc_vec[end] == curr_cyc) && (id_vec[end] == curr_id)) {
            end <- end + 1L
        }

        # the current cycle index are the indices between `start` and one before
        # `end`.  Updatate `start` and `ctr`.
        cyc_idx_list[[ctr]] <- start : (end - 1L)
        start <- end
        ctr <- ctr + 1L
    }

    # return the data shortened to the number of cycles
    cyc_idx_list[1:(ctr - 1L)]
}




# returns a vector of indices that maps the t-th element of `id` to the i-th unique value
#
# PRE: `id` is an atomic vector with length >= 1 that has already been sorted.

get_id_map <- function(id) {

    ctr <- 1L
    id_idx <- vector("integer", length(id))
    id_idx[1L] <- 1L

    for (i in seq_along(id_idx)[-1L]) {

        if (id[i] != id[i - 1]) {
            ctr <- ctr + 1L
        }

        id_idx[i] <- ctr
    }

    id_idx
}



# TODO: redo this back to -1 and 2 scoring

# TODO: I don't think these docs are current?

# converts the `daily$sex_yester` data to an integer vector of the same length
# as the data, where the elements take values -2L, -1L, 0L, and 1L, and which
# correspond to
#
#     -2L:  a missing value for sex yesterday where yesterday was the day before
#           the fertile window
#     -1L:  a missing value for sex yesterday where yesterday was one of the
#           days during the fertile window
#      0L:  sex yesterday in known to be no
#      1L:  sex yesterday is known to be yes
#
# PRE: daily is a data.frame with a column named "sex_yester" and a column for
# cycle day with name as given by `var_nm$fw`.  `fw_incl` is a vector with the
# values of the fertile window days, in the order that they occur.

get_sex_yester_coding <- function(daily, var_nm) {

    # case: no "sex yesterday" data in `daily`, i.e. we are not imputing missing
    # intercourse data
    if (! ("sex_yester" %in% colnames(daily))) {
        return(NULL)
    }

    # # hard-coded values which indicate what type of missingness occured
    # SEX_NO_INPUTE_BEFORE <- -2L
    # SEX_NO_INPUTE_DURING <- -1L

    # # map "sex yesterday" variable to binary data.  Missing values are
    # # preserved.
    # sex_yester <- map_vec_to_bool(daily$sex_yester) %>% as.integer
    # sex_yester_miss_idx <- which(is.na(sex_yester))

    # # bind `cycleday` variable to the fertile window day variable and set
    # # `fw_day_1` to be the value of the first day in the fertile window
    # cycleday <- daily[[var_nm$fw]]
    # fw_day_1 <- fw_incl[1L]

    # # each iteration checks one of the missing values of the "sex yesterday"
    # # variable and codes it so that the value of the variable indicates whether
    # # it occured on the day before the fertile window, or during the fertile
    # # window
    # for (i in sex_yester_miss_idx) {

    #     if (cycleday[i] == fw_day_1) {
    #         sex_yester[i] <- SEX_NO_INPUTE_BEFORE
    #     } else {
    #         sex_yester[i] <- SEX_NO_INPUTE_DURING
    #     }
    # }

    # sex_yester

    # map "sex yesterday" variable to binary data and where missing values are
    # coded as a -2, which is the code for an imputed value of "no sex
    # yesterday"
    sex_yester <- map_vec_to_bool(daily$sex_yester) %>% as.integer
    sex_yester[is.na(sex_yester)] <- -2L
    sex_yester
}




map_xmiss_to_x <- function(comb_dat, var_nm, xmiss) {

    # fills in the missing values of `xmiss` with the indices of the missing
    # intercourse.  This works because the same missing intercourse observations
    # are guaranteed to be in both vectors (and no more and no less), and the
    # order has been preserved.
    miss_bool <- is.na(xmiss)
    xmiss[miss_bool] <- which(is.na(comb_dat[[var_nm$sex]]))

    # subtract 1L for 0-based indexing
    xmiss - 1L
}
