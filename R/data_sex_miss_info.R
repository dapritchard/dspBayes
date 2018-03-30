get_sex_miss_info <- function(comb_dat, var_nm,  use_na) {

    # exit early if we're not imputing missing intercourse
    if (! ((use_na == "sex") || (use_na == "all"))) {
        list(x_map  = integer(0L),
             w_map  = integer(0L),
             xi_map = integer(0L))
    }

    # list with each element a vector of indices constituting a cycle
    cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)

    # whether this is a day during a pregnancy cycle
    preg_bool <- map_vec_to_bool(comb_dat[[var_nm$preg]])
    # whether intercourse was missing
    sex_miss_bool <- is.na(comb_dat[[var_nm$sex]])

    # a map from the daily data to pregnancy
    w_days_idx <- vector("integer", length(preg_bool))
    w_days_idx[preg_bool] <- seq_len(sum(preg_bool))

    # containers for the data we are constructing
    x_map  <- vector("integer", length(cyc_idx_list))
    w_map  <- vector("integer", length(cyc_idx_list))
    xi_map <- vector("integer", length(cyc_idx_list))

    # each iteration checks whether there is any missing intercourse in the
    # current cycle, and if so adds an entry to each of the mapping vectors
    ctr <- 1L
    for (curr_cyc_idx in cyc_idx_list) {

        # case: there is missing intercourse in this cycle, add an entry to our
        # mapping vectors
        if (any(sex_miss_bool[curr_cyc_idx])) {

            x_map[ctr]  <- curr_cyc_idx[1L]
            w_map[ctr]  <- w_days_idx[curr_cyc_idx[1L]]
            xi_map[ctr] <- xi_days_idx[curr_cyc_idx[1L]]

            ctr <- ctr + 1L
        }
    }

    # strip the unused portion of the vectors, and subtract by 1 for 0-based
    # indexing
    list(x_map  = x_map[seq_len(ctr - 1L)] - 1L,
         w_map  = w_map[seq_len(ctr - 1L)] - 1L,
         xi_map = xi_map[seq_len(ctr - 1L)] - 1L)
}
