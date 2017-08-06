# construct data ---------------------------------------------------------------

var_nm <- list(sex = "sex")
comb_dat <- data.frame(sex = c(1, 0, 1, NA, NA, NA),
                       age = c("35-37", "35-37", "<35", "<35", "38+", "38+"))

all_nonsex <- data.frame(sex = FALSE,
                         age = "35-37")


# begin testing ----------------------------------------------------------------

test_remove_days_no_sex <- function() {

    out <- remove_days_no_sex(comb_dat, var_nm)
    checkIdentical(comb_dat[-2, ], out)

    checkException(remove_days_no_sex(all_nonsex, var_nm))
}
