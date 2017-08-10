# construct data ---------------------------------------------------------------

daily <- data.frame(id  = rep(c(1, 2), each = 6),
                    cyc = rep(c(1, 2, 1, 2), each = 3),
                    day = rep(c(4, 2, 3), 4))

# arbitrary permutations of the data
perm_idx <- c(2, 10, 3, 1, 4, 12, 8, 9, 11, 5, 7, 6)
daily_perm <- daily[perm_idx, ]

# keep all variables
fw_incl <- c(4, 2, 3)

# `var_nm ` mock
var_nm <- list(id = "id", cyc = "cyc", fw = "day")




# begin testing ----------------------------------------------------------------

test_sort <- function() {

    # should be able to recover data before permutation
    checkEquals(daily, sort_dsp(daily, var_nm, fw_incl))

    # we should get the same result if the fertile window is days 2 and 3, and
    # the day before the fertile window is day 4
    checkEquals(daily, sort_dsp(daily, var_nm, c(2, 3), 4))
}
