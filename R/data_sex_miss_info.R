get_sex_miss_info <- function(comb_dat, var_nm,  use_na) {

    # exit early if we're not imputing missing intercourse
    if (! ((use_na == "sex") || (use_na == "all"))) {
        return(list(x_cyc    = integer(0L),
                    w_cyc    = integer(0L),
                    xi_idx   = integer(0L),
                    sex_prev = integer(0L),
                    n_cyc    = 0L))
    }

    # list with each element a vector of indices constituting a cycle
    cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)

    # whether this is a day during a pregnancy cycle
    preg_bool <- map_vec_to_bool(comb_dat[[var_nm$preg]])
    # whether intercourse was missing
    sex_miss_bool <- is.na(comb_dat[[var_nm$sex]])
    # whether intercourse occurred the day before.  Missing is coded as -99L.
    sex_yest_bin <- map_vec_to_bool(comb_dat$sex_yester) %>% as.integer
    is.na(sex_yest_bin) <- -99L

    # a map from the daily data to pregnancy
    w_days_idx <- vector("integer", length(preg_bool))
    w_days_idx[preg_bool] <- seq_len(sum(preg_bool))

    # a map from the daily data to subject
    id <- unique(comb_dat[[var_nm$id]])
    xi_days_idx <- rep(id, table(id))

    # containers for the data we are constructing
    x_cyc    <- vector("integer", length(cyc_idx_list))
    w_cyc    <- vector("integer", length(cyc_idx_list))
    xi_cyc   <- vector("integer", length(cyc_idx_list))
    sex_prev <- vector("integer", length(cyc_idx_list))

    # each iteration checks whether there is any missing intercourse in the
    # current cycle, and if so adds an entry to each of the mapping vectors
    ctr <- 1L
    for (curr_cyc_idx in cyc_idx_list) {

        # case: there is missing intercourse in this cycle, add an entry to our
        # mapping vectors
        if (any(sex_miss_bool[curr_cyc_idx])) {

            x_cyc[ctr]    <- curr_cyc_idx[1L]
            w_cyc[ctr]    <- w_days_idx[curr_cyc_idx[1L]]
            xi_cyc[ctr]   <- xi_days_idx[curr_cyc_idx[1L]]
            sex_prev[ctr] <- sex_yest_bin[curr_cyc_idx[1L]]

            ctr <- ctr + 1L
        }
    }

    # strip the unused portion of the vectors, and subtract by 1 for 0-based
    # indexing.  `n_cyc` given the length of the vectors (i.e. the number of
    # cycles)
    list(x_cyc    = x_cyc[seq_len(ctr - 1L)] - 1L,
         w_cyc    = w_cyc[seq_len(ctr - 1L)] - 1L,
         xi_idx   = xi_cyc[seq_len(ctr - 1L)] - 1L,
         sex_prev = sex_prev[seq_len(ctr - 1L)],
         n_cyc    = ctr - 1L)
}
