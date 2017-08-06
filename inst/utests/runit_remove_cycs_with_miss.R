# construct data ---------------------------------------------------------------

var_nm <- list(sex = "sex")
fw_incl <- 1:2

# corresponds to 3 cycles since the length of the fertile window is 2 days
sex <- c("y", "n", "y", NA, NA, NA)
age <- c(NA, NA, "35-37", "35-37", "38+", "38+")
comb_dat <- data.frame(sex = sex, age = age)

remove_cycs_with_miss_wrap <- function(use_na) {
    remove_cycs_with_miss(comb_dat, var_nm, fw_incl, use_na)
}


# begin testing ----------------------------------------------------------------

test_remove_cycs_with_miss <- function() {

    out <- remove_cycs_with_miss_wrap("all")
    checkIdentical(comb_dat, out)

    out <- remove_cycs_with_miss_wrap("intercourse")
    checkIdentical(comb_dat[3:6, ], out)

    out <- remove_cycs_with_miss_wrap("covariates")
    checkIdentical(comb_dat[1:2, ], out)

    checkException(remove_cycs_with_miss_wrap("none"))
}
