# construct data ---------------------------------------------------------------

base_df <- data.frame(matrix(1:4, 2, 2, dimnames = list(NULL, c("one", "two"))))
cyc_df <- data.frame(matrix(1:4, 2, 2, dimnames = list(NULL, c("a", "b"))))
day_df <- data.frame(matrix(1:4, 2, 2, dimnames = list(NULL, c("alpha", "beta"))))

no_nm <- base_df
colnames(no_nm) <- NULL


# begin testing ----------------------------------------------------------------

# check that we store the names of the special variables correctly
test_special_vars_nm <- function() {

    out <- consolidate_var_nm(formula(), NULL, NULL, NULL, "my_id", "", "", "")
    checkIdentical("my_id", out$id)

    out <- consolidate_var_nm(formula(), NULL, NULL, NULL, "", "my_cyc", "", "")
    checkIdentical("my_cyc", out$cyc)

    out <- consolidate_var_nm(formula(), NULL, NULL, NULL, "", "", "my_sex", "")
    checkIdentical("my_sex", out$sex)

    out <- consolidate_var_nm(formula(), NULL, NULL, NULL, "", "", "", "my_fw")
    checkIdentical("my_fw", out$fw)
}


test_extract_preg_nm <- function() {
    # the pregnancy variable is the outcome variable so it is extracted from the
    # LHS term in the model formula
    out <- consolidate_var_nm(alpha ~ beta, NULL, NULL, NULL, "", "", "", "")
    checkIdentical("alpha", out$preg)
}


test_extract_dataset_vars_nm <- function() {

    # if input is NULL for `baseline` or `cycle` then NULL is returned.  Note:
    # it is presumed that `daily` cannot be NULL.
    out <- consolidate_var_nm(alpha ~ beta + one, NULL, NULL, NULL, "w", "x", "y", "z")
    checkIdentical(character(0), out$all_base)
    checkIdentical(character(0), out$all_cyc)

    # we get `character(0)` if the data frames don't have column names
    out <- consolidate_var_nm(alpha ~ beta + one, no_nm, no_nm, no_nm, "w", "x", "y", "z")
    checkIdentical(character(0), out$all_base)
    checkIdentical(character(0), out$all_cyc)
    checkIdentical(character(0), out$all_day)

    # we get `character(0)` if none of the variables match
    out <- consolidate_var_nm(ten ~ eleven, base_df, cyc_df, day_df, "w", "x", "y", "z")
    checkIdentical(character(0), out$all_base)
    checkIdentical(character(0), out$all_cyc)
    checkIdentical(character(0), out$all_day)

    # check values for `baseline`
    out <- consolidate_var_nm(alpha ~ beta + two, base_df, NULL, NULL, "one", "", "", "")
    checkIdentical(c("one", "two"), out$all_base)

    # check values for `cycle`
    out <- consolidate_var_nm(alpha ~ beta + two, NULL, cyc_df, NULL, "a", "", "", "")
    checkIdentical("a", out$all_cyc)

    # check values for `daily`
    out <- consolidate_var_nm(alpha ~ beta + two, NULL, NULL, day_df, "", "", "", "")
    checkIdentical(c("alpha", "beta"), out$all_day)
}


test_extract_all_nm <- function() {

    # should have all four variables from the formula plus the four special
    # variable names, in some order
    out <- consolidate_var_nm(alpha ~ beta + one + b, NULL, NULL, NULL, "w", "x", "y", "z")
    checkTrue(all(c("alpha", "beta", "one", "b", "w", "x", "y", "z") %in% out$all))
    checkIdentical(8L, length(out$all))

    # if variable names are used twice then we consolidate them
    out <- consolidate_var_nm(alpha ~ beta + one + b, NULL, NULL, NULL, "one", "b", "y", "z")
    checkTrue(all(c("alpha", "beta", "one", "b", "y", "z") %in% out$all))
    checkIdentical(6L, length(out$all))
}
