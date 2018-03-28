# TODO: fcn description

get_xmiss <- function(comb_dat, var_nm, fw_incl, use_na) {

    # nothing to do if we're not imputing missing intercourse data
    if (! ((use_na == "sex") || (use_na == "all"))) {
        return(integer(0L))
    }

    # cycle length, number of cycles,
    K <- length(fw_incl)
    N <- K + 1L
    ncycs <- nrow(comb_dat) / K

    # indices in `comb_dat` for each cycle
    cyc_num_vec <- seq_len(ncycs) - 1L
    cyc_idx_list <- lapply(cyc_num_vec, function(i) ((i * K) + 1) : ((i + 1L) * K))

    # allocate the maximum possible amount of needed memory.  Each cycle with
    # missing intercourse is given 1 more day than the size of the fertile
    # window (the extra day is for intercourse the day before the start of the
    # fertile window).
    out_xmiss <- integer(ncycs * N)

    # standardize allowable forms of intercourse variable to -1/0 data.  Missing
    # values are mapped to the sentinal value of -98L.
    x_yest <- map_vec_to_bool(comb_dat$sex_yester) %>% as.integer %>% `-`(., 1L)
    x_yest[is.na(x_yest)] <- -98L

    # store intercourse data as a binary variable.  Missingness is preserved.
    X <- map_vec_to_bool(comb_dat[, var_nm$sex]) %>% as.integer %>% `-`(., 1L)

    # step through the cycles in the data and append the day before plus the FW
    # intercourse data if any intercourse is missing in the cycle.  Loop
    # invariant: `r` always points to the start of the unused memory in
    # `out_xmiss`.
    r <- 1L
    for (curr_idx in cyc_idx_list) {

        # intercourse values for current cycle
        day_before <- x_yest[curr_idx[1L]]
        curr_x <- X[curr_idx]

        # include current cycle if any intercourse data is missing
        if (any(is.na(curr_x))) {

            # copy day before and FW intercourse data to `out_xmiss`
            out_xmiss[r] <- day_before
            out_xmiss[(r + 1L) : (r + K)] <- curr_x

            # update loop invariant
            r <- r + N
        }
    }

    # subset to the size of the stored data
    out_xmiss[seq_len(r - 1L)]
}
