utest_xi <- function(dsp_data, W, ubeta, phi_val, seed_val) {

    # container to store samples in
    target_samples <- vector("numeric", length(dsp_data$subj_day_blocks))

    # create a vector indexing which subjects had a cycle with a pregnancy.  Add
    # 1 to subject indices since data is using 0-based indexing.
    w_subj_idx <- sapply(dsp_data$w_day_blocks, function(x) x["subj_idx"] + 1L)

    # tracks which index of W we are on
    w_ctr <- 1L

    # each iteration samples a value of xi for one subject
    set.seed(seed_val)
    for (i in seq_along(dsp_data$subj_day_blocks)) {

        # calculate `sum_jk W_ijk ` for current i
        curr_preg_cyc_idx <- which(i == w_subj_idx)
        if (length(curr_preg_cyc_idx) > 0L) {

            # find the right indices of W days and take the sum
            curr_w_day_block <- dsp_data$w_day_blocks[[curr_preg_cyc_idx]]
            curr_w_start_day <- w_ctr
            curr_w_end_day <- w_ctr + curr_w_day_block["n_days"] - 1L
            curr_w_sum <- W[curr_w_start_day : curr_w_end_day] %>% sum

            # point to the start of the next block of W values
            w_ctr <- w_ctr + curr_w_day_block["n_days"]
        }
        else {
            curr_w_sum <- 0
        }

        # calculate `sum_jk X_ijk * exp( u_{ijk}^T beta )`
        # add 1 to day indices since data is using 0-based indexing
        curr_subj_day_block <- dsp_data$subj_day_blocks[[i]]
        curr_start_day <- curr_subj_day_block["beg_idx"] + 1L
        curr_end_day <- curr_subj_day_block["beg_idx"] + curr_subj_day_block["n_days"]
        curr_day_idx <- curr_start_day : curr_end_day
        curr_exp_beta_sum <- sum(dsp_data$intercourse$X[curr_day_idx] * exp(ubeta[curr_day_idx]))

        # sample i-th value of xi
        target_samples[i] <- rgamma(1L, phi_val + curr_w_sum, phi_val + curr_exp_beta_sum)
    }

    target_samples
}
