# construct data ---------------------------------------------------------------

sex_numer <- c(1, 0, 1, NA, NA, NA)
sex_char <- c("y", "n", "y", NA, NA, NA)
sex_fact <- factor(c("y", "n", "y", NA, NA, NA))
sex_bool <- c(TRUE, FALSE, TRUE, NA, NA, NA)

target <- c(TRUE, FALSE, TRUE, NA, NA, NA)


# begin testing ----------------------------------------------------------------

map_vec_to_bool <- function() {

    out_numer <- map_vec_to_bool(sex_numer)
    checkIdentical(target, out_numer)

    out_char <- map_vec_to_bool(sex_char)
    checkIdentical(target, out_char)

    out_fact <- map_vec_to_bool(sex_fact)
    checkIdentical(target, out_fact)

    out_bool <- map_vec_to_bool(sex_bool)
    checkIdentical(target, out_bool)
}
