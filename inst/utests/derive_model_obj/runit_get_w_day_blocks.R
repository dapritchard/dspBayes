# construct data ---------------------------------------------------------------

preg <- c("no", "yes", "yes", "no", "no", "no", "yes", "yes", "yes", "yes")
cyc <- rep(c("first", "second", "fifth", "sixth"), c(1, 2, 3, 4))
id <- rep(c("abby", "jane"), c(3, 7))
comb_dat <- data.frame(id = id, cyc = cyc, preg = preg)

var_nm <- list(id = "id", cyc = "cyc", preg = "preg")

# use 0-based indexing for these values (not `n_days`, since it isn't an index)
target_preg_1 <- c(beg_idx = 1L, n_days = 2L, subj_idx = 0L, cyc_idx = 1L)
target_preg_2 <- c(beg_idx = 6L, n_days = 4L, subj_idx = 1L, cyc_idx = 3L)

# begin testing ----------------------------------------------------------------

test_get_w_day_blocks <- function() {

    out <- get_w_day_blocks(comb_dat, var_nm)

    checkIdentical(target_preg_1, out[[1L]])
    checkIdentical(target_preg_2, out[[2L]])
}
