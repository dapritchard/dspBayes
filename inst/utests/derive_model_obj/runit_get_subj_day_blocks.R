# construct data ---------------------------------------------------------------

comb_dat <- list(id = c("mary", "mary", "had", "a", "a", "a", "little", "lamb", "lamb"))
var_nm <- list(id = "id")

target <- list(c(beg_idx = 0L, n_days = 2L),
               c(beg_idx = 2L, n_days = 1L),
               c(beg_idx = 3L, n_days = 3L),
               c(beg_idx = 6L, n_days = 1L),
               c(beg_idx = 7L, n_days = 2L))


# begin testing ----------------------------------------------------------------

test_get_subj_day_blocks <- function() {

    out <- get_subj_day_blocks(comb_dat, var_nm)
    checkIdentical(target, out)
}
