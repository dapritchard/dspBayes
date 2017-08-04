# construct data ---------------------------------------------------------------

daily <- data.frame(id  = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4),
                    cyc = c(3, 3, 4, 4, 4, 4, 4, 1, 1, 2, 2, 4, 4, 4),
                    day = c(7, 5, 1, 4, 9, 6, 8, 4, 6, 5, 8, 1, 5, 8),
                    row = 1:14)

keypairs <- data.frame(id   = c(1, 1, 2, 2, 3, 4),
                       cyc  = c(3, 4, 1, 2, 1, 4))

fw_incl <- 5:8

vars_nm <- list(id = "id", cyc = "cyc", fw = "day")

min_days_req <- 0L




# begin testing ----------------------------------------------------------------

target <- data.frame(id  = rep(keypairs$id, each = 4),
                     cyc = rep(keypairs$cyc, each = 4),
                     day = rep(fw_incl, NROW(keypairs)),
                     row = c(2, NA, 1, NA, NA, 6, NA, 7, NA, 9, NA, NA, 10, NA,
                             NA, 11, NA, NA, NA, NA, 13, NA, NA, 14))

test_day_obs <- function() {

    # full data frame
    actual <- construct_day_obs(daily, keypairs, fw_incl, vars_nm, min_days_req)
    checkEquals(actual, target)

    # length-1 fertile window.  Set the fertile window to day 8.
    actual <- construct_day_obs(daily, keypairs, 8, vars_nm, min_days_req)
    target_length1 <- target[1:6*4, ]
    checkEquals(actual, target_length1, check.attributes = FALSE)

    # require 2 nonmissing days in the fertile window
    actual <- construct_day_obs(daily, keypairs, fw_incl, vars_nm, 2)
    target_length1 <- target[-c(9:12, 17:20), ]
    checkEquals(actual, target_length1, check.attributes = FALSE)
}


test_no_cycles <- function() {

    # require more nonmissing days than are in the fertile window
    msg <- "no cycles passed the criteria"
    checkException(construct_day_obs(daily, keypairs, fw_incl, vars_nm, 99), msg)
}
construct_day_obs(daily, keypairs, integer(0), vars_nm, 0)
