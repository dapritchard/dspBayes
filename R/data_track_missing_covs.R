get_cov_col_miss_info <- function(expanded_df, dsp_model, use_na) {

    if ((use_na != "covariates") && (use_na != "all")) {
        return(vector("list", 0L))
    }

    # the assign attribute of `model.matrix` is an integer vector with an entry
    # for each column in the design matrix, and corresponding to the term in the
    # formula (stored in `explan_var_nm`) which gave rise to the column.  A
    # value of 0 corresponds to the intercept, if any.
    map_expand_to_orig <- attr(expanded_df, "assign")
    unique_map_vals <- setdiff(map_expand_to_orig, 0L)

    # `terms()` provides a representation of the symbolic model that we can use
    # to access some of the information about the model
    dsp_terms <- terms(dsp_model)
    # variable names for the model RHS
    explan_var_nm <-  attr(dsp_terms, "term.labels")
    # matrix showing which variable in the original data that each variable in
    # the model is comprised of (i.e. shows where interaction terms come from,
    # etc)
    var_composition_mat <- attr(dsp_terms, "factors")
    var_orig_nm <- row.names(var_composition_mat)

    # categorical data variable names
    categ_nm <- attr(expanded_df, "contrasts") %>% names

    # container in which to place missing covariate information
    cov_miss_list <- vector("list", length(unique_map_vals))
    names(cov_miss_list) <- explan_var_nm

    # each iteration stores information for the k-th variable in a list
    # providing the indices in the design matrix that the variable corresponds
    # to, whether the variable is categorical, and if the variable has a
    # reference cell coding
    for (k in unique_map_vals) {

        cov_miss_list[[k]] <- vector("list", 0L)

        # current variable name and corresponding column indices in the design
        # matrix
        curr_var_nm <- explan_var_nm[k]
        curr_idx <- which(map_expand_to_orig == k)
        curr_data <- expanded_df[, curr_idx, drop = FALSE]
        cov_miss_list[[k]]$idx <- curr_idx

        # the names of the variables in the original data from which the current
        # design matrix variable is comprised
        curr_var_composition_nm <- var_orig_nm[ as.logical(var_composition_mat[, k]) ]
        cov_miss_list[[k]]$composition_nm <- curr_var_composition_nm

        # # store current variable column indices using 0-based indexing
        # cov_miss_list[[k]]$col_start <- head(curr_idx, 1L) - 1L
        # cov_miss_list[[k]]$col_end <- tail(curr_idx, 1L)

        # case: not a categorical variable
        if (! all(curr_var_composition_nm %in% categ_nm)) {
            cov_miss_list[[k]]$categ <- FALSE
            cov_miss_list[[k]]$all_cells <- TRUE
        }
        # case: a categorical variable.  Find out if it is also a reference cell
        # coding or not.
        else {

            cov_miss_list[[k]]$categ <- TRUE

            if (all(rowSums(curr_data) == 1L)) {
                cov_miss_list[[k]]$all_cells <- TRUE
            } else {
                cov_miss_list[[k]]$all_cells <- FALSE
            }
        }

        curr_n_categs <- length(curr_idx) + !cov_miss_list[[k]]$all_cells

        # calculate the empirical class distributions
        curr_empirical_probs <- vector("numeric", curr_n_categs)
        for (j in 1:(curr_n_categs - 1L)) {
            curr_empirical_probs[j] <- mean(expanded_df[, curr_idx[j]], na.rm = TRUE)
        }
        curr_empirical_probs[curr_n_categs] <- 1 - sum(curr_empirical_probs)
        cov_miss_list[[k]]$empirical_probs <- curr_empirical_probs

        cov_miss_list[[k]]$any_miss <- any(is.na(expanded_df[, curr_idx]))
    }

    Filter(function(x) x["any_miss"], cov_miss_list)
}




get_cov_row_miss_info <- function(comb_dat, var_nm, U, cov_col_miss_info) {

    # return without doing any work if there is no missing data
    if (length(cov_col_miss_info) == 0L) {
        return(vector("list", 0L))
    }

    # TODO: let's use and store these in the model object creation section and
    # pass them into the funciton
    subj_idx_list <- get_subj_idx_list(comb_dat, var_nm)
    cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)
    day_idx_list <- lapply(1:NROW(comb_dat), identity)

    # create a vector that maps days in the daily data to the corresponding
    # index the data that includes only days that occured in pregnancy cycles
    # (i.e. `W` days).  Days in the daily data that did not occur during a
    # pregnancy cycle are given a value of 0.  Note that we can assume that
    # there are no missing in the pregnancy variable because such data would
    # already have been removed.
    preg_idx <- vector("integer", NROW(comb_dat))
    preg_bool <- map_vec_to_bool( comb_dat[[var_nm$preg]] )
    preg_idx[preg_bool] <- seq_len( sum(preg_bool) )

    # similar to `preg_idx` but now with missing sex
    sex_miss_idx <- vector("integer", NROW(comb_dat))
    sex_miss_bool <- is.na( comb_dat[[var_nm$sex]] )
    sex_miss_idx[sex_miss_bool] <- seq_len( sum(sex_miss_bool) )

    # create a map to the missing covariate row indices.  The k-th element of
    # the `map_to_miss_idx` is 0 if there was no missing for any of the
    # covariates for this day.  If there was missing on the day then it provides
    # the number of days with missing covariates so far in the data.
    row_miss_bool <- !complete.cases(U)
    map_to_miss_idx <- vector("integer", NROW(comb_dat))
    map_to_miss_idx[row_miss_bool] <- seq_len( sum(row_miss_bool) )

    # similar to `map_to_miss_idx` but now for days in which there are both
    # missing covariates and missing intercourse status
    map_to_sex_idx <- vector("integer", NROW(comb_dat))
    cov_and_sex_miss_bool <- row_miss_bool & sex_miss_bool
    map_to_sex_idx[cov_and_sex_miss_bool] <- seq_len( sum(cov_and_sex_miss_bool) )

    # each iteration stores a list in the `i`-th element of `cov_row_miss_list`
    # providing the missing covariate information for the corresponding variable
    cov_row_miss_list <- vector("list", length(cov_col_miss_info))
    for (i in seq_along(cov_col_miss_info)) {

        curr_col_info <- cov_col_miss_info[[i]]

        curr_nm <- curr_col_info$composition_nm
        if (length(curr_nm) > 1L) {
            stop("the current version of this program does not support ",
                 "interaction terms that contain missing data", call. = FALSE)
        }

        # the blocks of data are determined by whether the variable is from the
        # baseline, cycle-specific, or day-specific data
        if (curr_nm %in% var_nm$all_base) {
            block_idx_list <- subj_idx_list
        } else if (curr_nm %in% var_nm$all_cyc) {
            block_idx_list <- cyc_idx_list
        } else {
            block_idx_list <- day_idx_list
        }

        cov_row_miss_list[[i]] <- get_miss_block_info(U,
                                                      curr_col_info,
                                                      block_idx_list,
                                                      map_to_miss_idx,
                                                      map_to_sex_idx,
                                                      sex_miss_bool,
                                                      get_id_map(comb_dat[[var_nm$id]]))
        # TODO: pass in `id_map` from elsewhere
    }

    list(cov_row_miss_list    = cov_row_miss_list,
         cov_miss_w_idx       = preg_idx[row_miss_bool],
         cov_miss_x_idx       = sex_miss_idx[row_miss_bool])
}





get_miss_block_info <- function(U,
                                cov_miss_info,
                                block_idx_list,
                                map_to_miss_idx,
                                map_to_sex_idx,
                                sex_miss_bool,
                                id_map) {

    # index in `U` of the first column in the block of columns corresponding the
    # the current unexpanded variable
    categ1_col_idx <- cov_miss_info$idx[1L]

    # container for the missing covariate block info
    miss_block_info <- vector("list", length(block_idx_list))
    ctr <- 1L

    # if there are any missing in the `curr_idx` then the iteration calculates
    # the necessary missing covariate information and stores it in the `ctr`-th
    # element of `miss_block_info`.  Otherwise, the iteration is a noop.
    for (curr_idx in block_idx_list) {

        # case: there are no missing values in this block of data, so continue
        # to the next block.  We need only check the first value since if any of
        # the elements in the block are missing, then all are missing.
        curr_first_idx <- curr_idx[1L]
        if (! is.na( U[curr_first_idx, categ1_col_idx] )) {
            next
        }

        # number of days in the current missing block, and the number of
        # days with missing intercourse in the current block
        curr_sex_bool <- sex_miss_bool[curr_idx]
        curr_n_days <- length(curr_idx)
        curr_n_sex_days <- sum(curr_sex_bool)

        # index to the set of data for days with both missing covariates and
        # missing intercourse status
        curr_sex_idx <- map_to_sex_idx[curr_idx]
        curr_beg_sex_idx <- ifelse(curr_n_sex_days > 0L,
                                   head(curr_sex_idx[curr_sex_bool], 1L),
                                   0L)

        # construct the missing block info needed for the sampler routines,
        # using 0-based indexing
        miss_block_info[[ctr]] <- c(beg_day_idx = curr_first_idx - 1L,
                                    n_days      = curr_n_days,
                                    beg_w_idx   = map_to_miss_idx[curr_first_idx] - 1L,
                                    beg_sex_idx = curr_beg_sex_idx - 1L,
                                    n_sex_days  = curr_n_sex_days,
                                    u_col       = categ1_col_idx - 1L,
                                    subj_idx    = id_map[curr_first_idx] - 1L)

        ctr <- ctr + 1L
    }

    miss_block_info[ seq_len(ctr - 1L) ]
}




get_ugen_info <- function(cov_col_miss_info, cov_row_miss_info) {

    ugen_info_list <- vector("list", length(cov_col_miss_info))

    for (i in seq_along(ugen_info_list)) {

        curr_col_info <- cov_col_miss_info[[i]]
        curr_row_miss_info <- cov_row_miss_info$cov_row_miss_list[[i]]

        max_n_days_miss <- sapply(curr_row_miss_info, function(x) x["n_days"]) %>% max
        max_n_sex_days_miss <- sapply(curr_row_miss_info, function(x) x["n_sex_days"]) %>% max

        # since we're switching to 0-based indexing, then
        # `tail(curr_col_info$idx, 1L)` gives us one past the last column index.
        # If however, we're using reference-cell coding, then we add another
        # column to represent the referernce category.
        col_start <- curr_col_info$idx[1L] - 1L
        col_end <- tail(curr_col_info$idx, 1L) + !curr_col_info$all_cells

        # if there is a reference category then we make it one past the last
        # column.  If not then just set the value to -1 as a flag.
        ref_col <- ifelse(curr_col_info$all_cells, -1, col_end)

        ugen_info_list[[i]] <- c(col_start           = col_start,
                                 col_end             = col_end,
                                 ref_col             = ref_col,
                                 n_categs            = col_end - col_start,
                                 max_n_days_miss     = max_n_days_miss,
                                 max_n_sex_days_miss = max_n_sex_days_miss,
                                 u_prior_probs       = cov_row_miss_info$empirical_probs[[i]],
                                 w_idx               = cov_row_miss_info$cov_miss_w_idx,
                                 x_idx               = cov_row_miss_info$cov_miss_x_idx)
    }

    ugen_info_list
}
