# construct data ---------------------------------------------------------------

# store as variables to ensure consistency between baseline, cycle, and daily
# datasets
id <- c(1, 2)
cyc <- c(1, 2, 1, 2)

# basic data variables
var_nm <- list(id = "id", cyc = "cyc", fw = "day", sex = "sex")

# construct daily data.  `daily` is without a pregnancy variable, and
# `daily_preg` has one added.
daily <- data.frame(id       = rep(id, each = 10),
                    cyc      = rep(cyc, each = 5),
                    day      = factor(rep(1:5, 4)),
                    sex      = rep(c("no", "yes", NA_character_), length.out = 20),
                    stringsAsFactors = FALSE)

# create a true `sex_yester` variable by shifting the data up one, and turning
# the days at the beginning of the fertile window into NAs
sex_yester <- daily$sex[c(1L, 1L:(NROW(daily) - 1L))]
sex_yester[c(1, 6, 11, 16] <- NA

# arbitrary permutations of the data
day_perm_idx <- c(14, 3, 17, 11, 6, 8, 1, 19, 5, 13, 2, 16, 7, 9, 10, 20, 15, 18, 4, 12)
daily_perm <- daily[day_perm_idx, ]


# begin testing ----------------------------------------------------------------

test_daily_add_sex_yester <- function() {

    out <- daily_add_sex_yester(daily, var_nm, c(2, 3, 4), 1)
    sorted_out <- out[order(out$id, out$cyc, out$day), ]

    fw_idx <- c(2:4, 7:9, 12:14, 17:19)
    non_fw_idx <- setdiff(1:20, fw_idx)

    target <- data.frame(daily, sex_yester, stringsAsFactors = FALSE)

    # the rows should be identical in the fertile window without sorting
    checkIdentical(target[fw_idx, ], out[fw_idx, ])

    # the output should be the same as original data up to row permutations, not
    # including the added column
    checkIdentical(daily, sorted_out[, ("sex_yester" != colnames(sorted_out))])

    # the terms that are outside of the fertile window should be NAs
    checkTrue(out$sex_yester[non_fw_idx] %>% is.na %>% all)
}
