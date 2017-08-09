# conditionally remove any cycles from `comb_dat` that have either missing data
# in the intercourse variable, other covariates, or both, depending on the value
# of `use_na`
#
# PRE: `comb_data` is a data frame with column names as specified by `var_nm`.
# `fw_incl` is a vector specifying the days in the fertile window, and `use_na`
# is a string specifying the type of missingness to allow.

remove_cycs_with_miss <- function(comb_dat, var_nm, fw_incl, use_na) {

    if (identical(use_na, "all")) {
        return(comb_dat)
    }

    # fertile window length and number of cycles in the data
    fw_len <- length(fw_incl)
    n_cyc <- NROW(comb_dat) / fw_len

    # indices of the variables which must be nonmissing for us to keep the cycle
    # case: use missing for intercourse, but not for covariates
    if (identical(use_na, "intercourse")) {
        check_idx <- which(colnames(comb_dat) != var_nm$sex)
    }
    # case: use missing for covariates, but not for intercourse
    else if (identical(use_na, "covariates")) {
        check_idx <- which(colnames(comb_dat) == var_nm$sex)
    }
    # case: don't use any cycles with any missing
    else if (identical(use_na, "none")) {
        check_idx <- seq_along(comb_dat)
    }

    # container to store the indices of the columns that we wish to keep,
    # depending on whether and what type of missingness is in the cycle
    keep_idx_list <- vector("list", n_cyc)

    # each iteration in the loop checks the (k + 1)-th cycle for any missing in
    # the columns indexed by `check_idx`, and if none are missing then stores
    # the indices corresponding to the cycle in `keep_idx_list`
    for (k in (seq_len(n_cyc) - 1L)) {

        # obtain data in current cycle
        curr_idx <- (k * fw_len + 1L) : ((k + 1L) * fw_len)
        curr_df <- comb_dat[curr_idx, check_idx, drop = FALSE]

        # if all of the observations in the cycle are nonmissing for the data
        # columns of interest, then store the indices corresponding to the cycle
        # in `keep_idx_list`
        if (complete.cases(curr_df) %>% all) {
            keep_idx_list[[k + 1L]] <- curr_idx
        }
    }

    # flatten all of the indices stored in `keep_idx_list` into an atomic
    # vector, and ensure that there are at least some cycles without any missing
    keep_idx <- unlist(keep_idx_list)
    if (is.null(keep_idx)) {
        stop("no cycles existed without any missing", call. = FALSE)
    }

    comb_dat[keep_idx, , drop = FALSE]
}




# takes an atomic vector of type logical, numerical, factor, or character and
# maps it to a logical vector.  Missing values are preserved.

map_vec_to_bool <- function(sex_vec) {

    allowed_yes <- c("y", "Y", "yes", "Yes", "YES")

    if (is.logical(sex_vec)) {
        return(sex_vec)
    }
    else if (is.numeric(sex_vec)) {
        return(sex_vec != 0)
    }
    else if (is.factor(sex_vec) || is.character(sex_vec)) {
        out <- rep(NA, length(sex_vec))
        out[! is.na(sex_vec)] <- sex_vec[! is.na(sex_vec)] %in% allowed_yes
        return(out)
    }
    else {
        # illegal form of intercourse status
        stop("invalid form of sex_vec")
    }
}
