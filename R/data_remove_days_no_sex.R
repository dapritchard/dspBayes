# removes all observations from `comb_dat` for which intercourse did not occur.
# Observations with a missing value for intercourse remain in the data.
#
# PRE: `comb_dat` is a data frame with one of the column names given by
# `var_nm$sex`.

remove_days_no_sex <- function(comb_dat, var_nm) {

    # obtain the column index for the intercourse variable
    sex_col_idx <- which(colnames(comb_dat) == var_nm$sex)

    # obtain a Boolean vector indicating the days in the data for which
    # intercourse occured.  Days in the data for which intercourse was missing
    # will have a missing value for this vector also.
    yes_sex_bool <- map_vec_to_bool(comb_dat[[sex_col_idx]])

    # keep all days in which intercourse occured or was missing
    keep_idx <- which(is.na(yes_sex_bool) | yes_sex_bool)
    if (identical(length(keep_idx), 0L)) {
        stop("there were no days in which intercourse occured", call. = FALSE)
    }

    # return data subset to only days in which intercourse occured or was
    # missing
    comb_dat[keep_idx, , drop = FALSE]
}
