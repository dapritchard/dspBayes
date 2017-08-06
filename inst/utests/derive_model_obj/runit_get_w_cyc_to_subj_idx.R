# construct data ---------------------------------------------------------------

w_day_blocks <- list(c(beg_idx =  1L, n_days = 2L, subj_idx = 0L, cyc_idx = 1L),
                     c(beg_idx =  6L, n_days = 4L, subj_idx = 1L, cyc_idx = 3L),
                     c(beg_idx = 18L, n_days = 3L, subj_idx = 3L, cyc_idx = 6L))

target <- c(0L, 1L, 3L)


# begin testing ----------------------------------------------------------------

test_get_subj_idx <- function() {

    out <- get_w_cyc_to_subj_idx(w_day_blocks)
    checkIdentical(target, out)
}
