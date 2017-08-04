merge_dsp_data <- function(baseline, cycle, daily, var_nm, fw_incl, min_days_req) {

    # reduce each dataset to only contain the variables that we need.  Returns
    # NULL if the input is NULL.
    base_red <- get_red_dat(baseline, var_nm)
    cyc_red <- get_red_dat(cycle, var_nm)
    day_red <- get_red_dat(daily, var_nm)

    # obtain every unique subject id and cycle pair for any cycle that has
    # pregnancy outcome data
    keypairs <- get_keypairs(day_red, cyc_red, var_nm)

    # reduce daily dataset to only include rows corresponding to days in the
    # fertile window.  Note that we must do this after obtaining `keypairs`,
    # because otherwise we could lose cycles that have pregnancy information but
    # where none of the fertile window days are recorded.
    keep_idx <- day_red[[var_nm$fw]] %in% fw_incl %>% which
    day_red <- day_red[keep_idx, ]

    # a data frame with observations given by the cross product of every day in
    # the fertile window and every (id, cycle) keypair
    comb_dat <- construct_day_obs(day_red, keypairs, fw_incl, var_nm, min_days_req)

    # join operation on id and cycle.  `all.x = TRUE` specifies that we keep
    # every observation in `day_formatted`.  Under these specifications `merge`
    # throws away observations in `cyc_red` that don't have matches, which is
    # what we want since `day_formatted` already has all of the id and cycle
    # pairs with pregnancy information.  We defer sorting the joined data until
    # later.
    if (! is.null(cyc_red)) {
        # TODO: check that there are not matches among model variable names
        # other than the ones being joined against
        comb_dat <- merge(x     = comb_dat,
                          y     = cyc_red,
                          by    = c(var_nm$id, var_nm$cyc),
                          all.x = TRUE,
                          sort  = FALSE)
    }

    # conditionally merge baseline in with cycle and daily.  See above merge for
    # discussion of the parameter settings.
    if (! is.null(base_red)) {
        # TODO: check that there are not matches among model variable names
        # other than the ones being joined against
        comb_dat <- merge(x     = comb_dat,
                          y     = base_red,
                          by    = var_nm$id,
                          all.x = TRUE,
                          sort  = FALSE)
    }

    # drop unused factor levels after subsetting
    for (i in seq_along(comb_dat)) {
        if (is.factor(comb_dat[[i]])) {
            comb_dat[[i]] <- droplevels(comb_dat[[i]])
        }
    }

    # return data after sorting it
    sort_dsp(comb_dat, keypairs, fw_incl, var_nm)
}




# sort over id / cycle / fertile window
sort_dsp <- function(comb_dat, keypairs, fw_incl, var_nm) {
    # turn fertile window vector into a factor if it's not already so that
    # `order` will give us the desired sorting
    out <- comb_dat[order(comb_dat[, var_nm$id],
                          comb_dat[, var_nm$cyc],
                          comb_dat[, var_nm$fw] %>% factor(., levels = fw_incl)), ]
    row.names(out) = NROW(comb_dat) %>% seq_len
    out
}




# takes a data frame `dataset` and returns the data after removing any columns
# with names do not appear in `var_nm$all`
#
# PRE: `dataset` is a data.frame and `var_nm` is a list with an element `all`
# which is a character vector

get_red_dat <- function(dataset, var_nm) {

    # if data is NULL then function is a noop
    if (is.null(dataset)) {
        return(NULL)
    }

    # subset the data variables according to whether they are among the set of
    # all variables listed in the model or used to combine the data
    var_incl_bool <- colnames(dataset) %in% var_nm$all

    # we can't allow missing data for these fundamental columns, so create an
    # index of rows to remove
    critical_cols_bool <- colnames(dataset) %in% c(var_nm$id, var_nm$cyc, var_nm$fw, var_nm$preg)
    obs_incl_bool <- complete.cases(dataset[, critical_cols_bool])

    dataset[obs_incl_bool, var_incl_bool, drop = FALSE] %>%
        as.data.frame(., stringsAsFactors = FALSE)
}




# return a data frame with two columns, with the rows containing every unique
# subject id and cycle pair for any cycle that has pregnancy outcome data.
# `cycle` is allowed to be NULL.
#
# PRE: assumes `daily` and `cycle` are data frames (or cycle may be NULL), each
# having columns for id and cycle.  Either one or both of `daily` or `cycle`
# must also have a column for pregnancy status.  The names of each of these
# columns is given by `var_nm`.

get_keypairs <- function(daily, cycle, var_nm) {

    # create `keypairs_df` with all of the subject id and cycle pair for any
    # cycle that has pregnancy outcome data

    # case: both datasets contain pregnancy information
    if ((var_nm$preg %in% colnames(daily)) && (var_nm$preg %in% colnames(cycle))) {
        keypairs_df <- rbind(cycle[, c(var_nm$id, var_nm$cyc)],
                             daily[, c(var_nm$id, var_nm$cyc)],
                             stringsAsFactors = FALSE)
    }
    # case: only the daily data constains pregnancy information
    else if (var_nm$preg %in% colnames(daily)) {
        keypairs_df <- daily[, c(var_nm$id, var_nm$cyc)]
    }
    # case: only the cycle data constains pregnancy information
    else if (var_nm$preg %in% colnames(cycle)) {
        keypairs_df <- cycle[, c(var_nm$id, var_nm$cyc)]
    }
    # case: can't find the pregnancy status variable
    else {
        paste0("cannot find ", var_nm$preg, " in cycle or daily") %>% stop(call. = FALSE)
    }

    # obtain unique keypairs
    dup_keys_bool <- duplicated(keypairs_df)
    keypairs_df <- keypairs_df[! dup_keys_bool, ]

    # return keypairs sorted by id and cycle
    keypairs_df[order(keypairs_df[, var_nm$id], keypairs_df[, var_nm$cyc]), ]
}




# construct the day-specific part of the data.  Returns a data frame with
# observations given by the cross product of every day in the fertile window and
# every (id, cycle) keypair.
#
# PRE: assumes `daily` is a data frame with id, cycle and fertile window columns
# with names as given in `var_nm`.  `keypairs` is a nonempty data frame with
# exactly two columns for id and cycle.  `fw_incl` is a nonempty atomic
# vector. `var_nm` is a list providing the names of the id, cycle, and fertile
# window columns.  `min_days_req` is a length-1 numeric vector.

construct_day_obs <- function(daily, keypairs, fw_incl, var_nm, min_days_req) {

    # the length of fertile window and the number of cycles in the data
    fw_len <- length(fw_incl)
    n_cycle <- NROW(keypairs)

    # a data frame with vectors of the same types as in `daily` and number of
    # rows the same as the number of fertile window days.  All values are NA.
    template_df <- construct_template_df(daily, fw_incl)

    # container with each element a data frame representing the days in the
    # fertile window
    days_by_cyc <- vector("list", n_cycle)
    n_obs_by_cyc <- vector("logical", n_cycle)

    # each iteration constructs a data frame for the daily data corresponding to
    # current subject and cycle, and saves results to `days_by_cyc`
    for (k in 1:n_cycle) {

        # subject id and cycle for current keypair
        curr_id <- keypairs[k, var_nm$id]
        curr_cyc <- keypairs[k, var_nm$cyc]

        # tracks how many observations are nonmissing for the current subject and
        # cycle
        curr_n_obs <- 0L

        # copy empty data frame and fill in id, cycle, and fertile window day
        # data.  Note that `curr_id` and `curr_cyc` are scalars and rely on R to
        # recycle their values to fill the vector, while `fw_incl` is of the
        # appropriate length.
        curr_df <- template_df
        curr_df[[var_nm$id]] <- curr_id
        curr_df[[var_nm$cyc]] <- curr_cyc
        curr_df[[var_nm$fw]] <- fw_incl

        # reduce to daily to the current subject and cycle
        curr_day_idx <- (daily[[var_nm$id]] == curr_id &
                         daily[[var_nm$cyc]] == curr_cyc) %>% which
        curr_daily <- daily[curr_day_idx, ]

        # each iteration copies one row of data from `curr_daily` into the
        # correct row of `curr_df` if one exists, otherwise nothing is done
        for (i in seq_along(fw_incl)) {

            # row in `curr_daily` corresponding the fertile window day
            row_idx <- which(curr_daily[[var_nm$fw]] == fw_incl[i])

            # case: found exactly 1 match, copy contents from daily data into
            # current data frame
            if (length(row_idx) == 1L) {
                curr_df[i, ] <- curr_daily[row_idx, ]
                curr_n_obs <- curr_n_obs + 1L
            }
            # case: multiple matches, throw an error
            else if (length(row_idx) > 1L) {
                paste0("multiple matches in daily data for subject ", curr_id,
                       " cycle ", curr_cyc, " fw day ", fw_incl[i],
                       ": dropping observation") %>% stop(call. = FALSE)
            }
            # else: 0 matches found.  Do nothing (i.e. leave an observation with
            # missing values in it)

        } # end copy day-specific data loop

        # save contents of daily data corresponding to current subject and cycle
        # to storage list
        days_by_cyc[[k]] <- curr_df
        n_obs_by_cyc[k] <- curr_n_obs

    } # end construct daily data for a given keypair loop

    # remove observations that don't meet the minimum number of observations
    # specified to keep
    days_by_cyc <- days_by_cyc[n_obs_by_cyc >= min_days_req]

    # rbind each of the elements in `days_by_cyc`, each of which are data frames
    # corresponding to the days for a given subject and cycle pair
    rbind_similar_dfs(days_by_cyc)
}




# construct a data frame with vectors of the same types and attributes as
# `daily`, number of rows the same length as `fw_incl`, and all values are NA
#
# PRE: `daily` is a data frame of at least one row and one column, and `fw_incl`
# is a nonempty atomic vector with length no bigger than the number of rows in
# `daily`

construct_template_df <- function(daily, fw_incl) {

    # obtain a data frame of the right size and vector types
    template_df <- daily[seq_along(fw_incl), , drop = FALSE]

    # fill all of the values with missing
    for (k in seq_along(template_df)) {
        # index all of the vector so that R keeps the vector internal structure
        # and attributes
        template_df[[k]][seq_along(fw_incl)] <- NA
    }

    row.names(template_df) <- seq_along(fw_incl)
    template_df
}




# combine a list of data frames that have a specific form.  It is assumed that
# each element of `list_of_dfs` is a data frame of exactly the same dimensions
# and attributes.  The function then performs an operation that is logically
# equivalent to `rbind()`ing each of these data frames
#
# PRE: `list_of_dfs` is a list with length > 0 such that each element is a data
# frame of exactly the same dimensions and attributes

rbind_similar_dfs <- function(list_of_dfs) {

    if (length(list_of_dfs) == 0L) {
        stop("no cycles passed the criteria", call. = FALSE)
    }

    # set data parameters
    n_dfs <- length(list_of_dfs)
    n_vars <- list_of_dfs[[1L]] %>% NCOL
    n_elem <- list_of_dfs[[1L]] %>% NROW
    vec_len <- n_dfs * n_elem

    # container to put each of the concatenated columns into
    combined_list <- vector("list", n_vars)
    names(combined_list) <- list_of_dfs[[1L]] %>% names

    # each iteration creates one variable that is the concatenation of the j-th
    # variable across all of the data frames
    for (j in seq_len(n_vars)) {

        # initialize current vector in the list
        combined_list[[j]] <- rep(list_of_dfs[[1L]][, j], n_dfs)

        # construct index { 2, 3..., n_dfs }.  Note that this results in an
        # integer(0) vector if `n_dfs` equals 1.
        df_idx <- seq_len(n_dfs)
        df_idx <- df_idx[-n_dfs] + 1L

        # each iteration fills in the values of the j-th variable of the k-th
        # data frame into the appropriate elements of the j-th variable in
        # `combined_list`
        for (k in df_idx) {

            # indices of elements to fill in for current vector
            curr_idx <- ((k - 1) * n_elem + 1) : (k * n_elem)

            # copy the j-th column of the k-th data frame into the appropriate
            # elements of the combined_list column
            combined_list[[j]][curr_idx] <- list_of_dfs[[k]][, j]
        }
    }

    data.frame(combined_list, stringsAsFactors = FALSE)
}
