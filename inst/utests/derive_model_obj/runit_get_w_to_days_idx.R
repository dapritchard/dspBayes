# construct data ---------------------------------------------------------------

comb_dat <- data.frame(preg = c("no", "yes", "yes", "no", "no", "no", "yes", "yes", "yes", "yes"))
var_nm <- list(preg = "preg")

# use 0-based indexing for these values
target <- c(1:2, 6:9, -1L)


# begin testing ----------------------------------------------------------------

test_get_w_days_idx <- function() {

    out <- get_w_to_days_idx(comb_dat, var_nm)
    checkIdentical(target, out)
}
