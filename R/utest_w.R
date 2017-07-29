utest_w <- function(dsp_data, xi, ubeta, seed_val) {

    # sample from a zero-truncated Poisson distribution
    rpois_zero_tr <- function(lambda) {
        runif(length(lambda), exp(-lambda), 1) %>% qpois(., lambda)
    }

    # for convenience
    w_day_blocks <- dsp_data$w_day_blocks
    X <- dsp_data$X

    # container to store samples in
    w_list <- vector("list", length(w_day_blocks))

    # each iteration stores a vector of W samples in the i-th element of
    # `w_list`.  Each element of the list is for one pregnancy cycle, and each
    # element in the vectors is one day in the cycle.
    set.seed(seed_val)
    for (i in seq_along(w_day_blocks)) {

        # find the indices in the day-specific data corresponding to the current
        # pregnancy cycle.  Add 1 to adjust for 0-based indexing.
        curr_block <- w_day_blocks[[i]]
        curr_idx <- (curr_block["beg_idx"] + 1L) : (curr_block["beg_idx"] + curr_block["n_days"])
        curr_subj <- curr_block["subj_idx"] + 1L

        # calculate the means for the multinomial and Poisson densities.  The
        # multinomial means are unnormalized (i.e. don't sum to 1).
        day_means <- X[curr_idx] * exp(ubeta[curr_idx])
        pois_mean <- xi[curr_subj] * sum(day_means)

        # sample `sum_k W_ijk`, and then sample `W_ij | sum_k W_ijk`
        w_sum <- rpois_zero_tr(pois_mean)
        w_list[[i]] <- rmultinom(1L, w_sum, day_means) %>% drop
    }

    unlist(w_list)
}
