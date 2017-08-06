# construct data ---------------------------------------------------------------

subj_day_blocks <- list(c(beg_idx = 0L, n_days = 2L),
                        c(beg_idx = 2L, n_days = 1L),
                        c(beg_idx = 3L, n_days = 3L),
                        c(beg_idx = 6L, n_days = 1L),
                        c(beg_idx = 7L, n_days = 2L))

target <- c(rep(0L, 2L), 1L, rep(2L, 3L), rep(3L, 1L), rep(4L, 2L))


# begin testing ----------------------------------------------------------------

test_get_day_to_subj_idx <- function() {

    out <- get_day_to_subj_idx(subj_day_blocks)
    checkIdentical(target, out)
}
