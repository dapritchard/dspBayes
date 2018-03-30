strip_days <- function(comb_dat, var_nm, fw_incl, use_na) {

    # early exit if we don't allow any missing
    if (use_na == "none") {
        return(list(comb_dat         = comb_dat,
                    miss_sex_cyc_idx = integer(0L),
                    sex_yest         = integer(0L)))
    }

    # cycle length, number of cycles,
    K <- length(fw_incl)
    ncycs <- nrow(comb_dat) / K

    # indices in `comb_dat` for each cycle
    cyc_num_vec <- seq_len(ncycs) - 1L
    cyc_idx_list <- lapply(cyc_num_vec, function(i) ((i * K) + 1) : ((i + 1L) * K))

    # indices to the first day of a cycle with missing intercourse and
    # intercourse status for the 0-th day of such a cycle
    sex_miss_bool <- vector("logical", nrow(comb_dat))
    miss_sex_cyc_idx <- vector("integer", length(cyc_num_vec))
    day_zero_sex <- vector("integer", length(cyc_num_vec))
    ctr <- 1L

    # each iteration checks whether the current cycle has any missing
    # intercourse days in it, and if so sets the corresponding days in
    # `sex_miss_bool` to TRUE.  Additionally, the index of the first day in the
    # cycle is recorded in `miss_sex_cyc_idx`, and the intercourse status for
    # the 0-th day of the cycle is recorded in `day_zero_sex`
    sex_yest_binary <- map_vec_to_bool(comb_dat$sex_yester) %>% as.integer
    for (curr_cyc_idx in cyc_idx_list) {

        curr_sex <- comb_dat[curr_cyc_idx, var_nm$sex]

        if (any(is.na(curr_sex))) {

            sex_miss_bool[curr_cyc_idx] <- TRUE
            miss_sex_cyc_idx[ctr] <- curr_cyc_idx[1L]
            day_zero_sex[ctr] <- sex_yest_binary[curr_cyc_idx[1L]]
            ctr <- ctr + 1L
        }
    }

    # strip unused portion of vectors
    miss_sex_cyc_idx <- miss_sex_cyc_idx[seq_len(ctr - 1L)]
    day_zero_sex <- day_zero_sex[seq_len(ctr - 1L)]

    # whether any covariates are missing
    miss_covs_bool <- !complete.cases(comb_dat[, var_nm$covs])

    # whether intercourse occured on the given day (or was missing, but this
    # will get picked up by `sex_miss_bool` anyway)
    sex_yes_bool <- map_vec_to_bool(comb_dat[[var_nm$sex]])
    sex_yes_bool[is.na(sex_yes_bool)] <- TRUE

    # data where any sex occured, or else we need to keep the day because there
    # were missing covariates or it was in a cycle with missing intercourse
    comb_dat <- comb_dat[sex_yes_bool | miss_covs_bool | sex_miss_bool, ]

    list(comb_dat         = comb_dat,
         miss_sex_cyc_idx = miss_sex_cyc_idx,
         day_zero_sex     = day_zero_sex)
}
