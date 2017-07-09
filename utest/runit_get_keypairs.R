# construct data ---------------------------------------------------------------

daily <- data.frame(id   = c(1, 1, 1, 1, 3, 3),
                    cyc  = c(1, 1, 2, 2, 1, 1),
                    preg = rep("no", 6))

cycle <- data.frame(id   = c(1, 1, 8, 8, 9, 9),
                    cyc  = c(1, 2, 1, 2, 1, 2),
                    preg = rep("no", 6))

vars_nm <- list(id = "id", cyc = "cyc", preg = "preg")

target_all <- data.frame(id  = c(1, 1, 8, 8, 9, 9, 3),
                         cyc = c(1, 2, 1, 2, 1, 2, 1))




# begin testing ----------------------------------------------------------------

test_valid_inp <- function() {

    # both daily and cycle
    actual <- get_keypairs(daily, cycle, vars_nm)
    checkEquals(target_all, actual, check.attributes = FALSE)

    # daily has the pregnancy information
    target <- target_all[c(1, 2, 7), ]
    actual <- get_keypairs(daily, cycle[, -3L], vars_nm)
    checkEquals(target, actual, , check.attributes = FALSE)

    # cycle is NULL
    target <- target_all[c(1, 2, 7), ]
    actual <- get_keypairs(daily, NULL, vars_nm)
    checkEquals(target, actual, , check.attributes = FALSE)

    # cycle has the pregnancy information
    target <- target_all[1:6, ]
    actual <- get_keypairs(daily[, -3L], cycle, vars_nm)
    checkEquals(target, actual, , check.attributes = FALSE)
}
