# construct data ---------------------------------------------------------------

A <- data.frame(matrix(1:4, 2, 2, dimnames = list(NULL, c("one", "two"))))
zero_col_df <- A[, integer(0)]
one_col_df <- A[, 1L, drop = FALSE]
two_col_df <- A[, 1:2]

var_nm_none_match <- list(all = character(0))
var_nm_one_match <- list(all = c("zero", "one", "three"))
var_nm_all_match <- list(all = c("zero", "one", "two", "three"))


# begin testing ----------------------------------------------------------------

test_inputs <- function() {

    # NULL input
    checkIdentical(NULL, get_red_dat(NULL, var_nm_all_match))

    # should get back the same data frame if all vars are in the model
    checkIdentical(zero_col_df, get_red_dat(zero_col_df, var_nm_all_match))
    checkIdentical(one_col_df, get_red_dat(one_col_df, var_nm_all_match))
    checkIdentical(two_col_df, get_red_dat(two_col_df, var_nm_all_match))

    # now match against 1 variable name
    checkIdentical(one_col_df, get_red_dat(one_col_df, var_nm_one_match))
    checkIdentical(one_col_df, get_red_dat(two_col_df, var_nm_one_match))

    # now match against 0 variables name
    checkIdentical(zero_col_df, get_red_dat(two_col_df, var_nm_none_match))
}
