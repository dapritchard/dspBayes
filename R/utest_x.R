utest_x <- function(dsp_data, W, xi, ubeta, utau, seed_val) {

    sample <- function(W, xi, ubeta, utau) {

        w_ctr <- 0L
        for (miss_cyc in x_miss_cyc_list) {

            sample_cycle(miss_cyc, W, xi, ubeta, utau)
        }
    }


    sample_cycle <- function(miss_cyc, W, xi, ubeta, utau) {

        miss_day <- x_miss_day_list[[miss_cyc["beg_idx"] + 1L]]
        xi_i <- xi[miss_cyc["subj_idx"] + 1L]

        if (miss_day["prev"] < 0L) {
            prev_day_sex <- sample_day_before_fw_sex()
        }

        w_idx <- miss_cyc["preg_idx"]
        cycle_miss_day_idx <- seq(miss_cyc["beg_idx"] + 1L,
                                  miss_cyc["beg_idx"] + miss_cyc["n_days"],
                                  1L)

        for (curr_miss_day_idx in cycle_miss_day_idx) {

            miss_day <- x_miss_day_list[[curr_miss_day_idx]]
            curr_day_idx <- miss_day["idx"] + 1L

            if (w_idx != -1L) {

                w_idx <- w_idx + 1L

                if (W[w_idx] > 0L) {
                    X[curr_day_idx] <<- prev_day_sex <- 1L
                    next
                }
            }

            if (miss_day["prev"] >= 0L) {
                prev_day_sex <- miss_day["prev"]
            }

            prior_prob_yes <- calc_prior_prob(utau, curr_miss_day_idx, prev_day_sex)
            posterior_prob_yes <- calc_posterior_prob(ubeta, xi_i, curr_day_idx)

            X[curr_day_idx] <<- prev_day_sex <- sample_x_ijk(prior_prob_yes,
                                                             posterior_prob_yes)
        }
    }


    calc_prior_prob <- function(utau, miss_day_idx, prev_day_sex) {
        if (prev_day_sex == 0L) {
            1 / (1 + exp(-utau[miss_day_idx]))
        } else {
            1 / (1 + exp(-utau[miss_day_idx] - sex_coef))
        }
    }


    calc_posterior_prob <- function(ubeta, xi_i, day_idx) {

        exp(-xi_i * exp(ubeta[day_idx]))
    }


    sample_x_ijk <- function(prior_prob_yes, posterior_prob_yes) {

        unnormalized_prob_no <- 1 - prior_prob_yes
        unnormalized_prob_yes <- prior_prob_yes * posterior_prob_yes
        prob_no <- unnormalized_prob_no / (unnormalized_prob_no + unnormalized_prob_yes)

        as.integer(runif(1) > prob_no)
    }


    sample_day_before_fw_sex <- function() {
        as.integer(runif(1L) < cohort_sex_prob)
    }


    # function constants.  Arbitrary choice for `N_TEST_SAMPLES`.
    SEX_MISS <- 2L
    N_TEST_SAMPLES <- 20L

    # bind variable names to the data for convenience
    X <- dsp_data$intercourse$X
    x_miss_cyc_list <- dsp_data$intercourse$miss_cyc
    x_miss_day_list <- dsp_data$intercourse$miss_day
    sex_coef <- dsp_data$tau_fit$sex_coef
    cohort_sex_prob <- dsp_data$tau_fit$cohort_sex_prob

    # arbitrary choices of missing day, prior / posterior probabilities, and
    # xi_i
    test_data_miss_day_idx <- 1L
    test_data_miss_day <- x_miss_day_list[[test_data_miss_day_idx]]
    test_data_nonmiss_day_idx <- unname(test_data_miss_day["idx"] + 1L)
    test_data_prior_prob <- 0.5
    test_data_posterior_prob <- 0.9
    test_data_xi_i <- 1.2

    # update X
    set.seed(seed_val)
    sample(W, xi, ubeta, utau)

    # calculate prior probabilities for `X_ijk`
    target_prior_prob_no_prev <- calc_prior_prob(utau, test_data_miss_day_idx, 0L)
    target_prior_prob_yes_prev <- calc_prior_prob(utau, test_data_miss_day_idx, 1L)

    # calculate posterior probabilities for `X_ijk`
    target_posterior_prob <- calc_posterior_prob(ubeta, test_data_xi_i, test_data_nonmiss_day_idx)

    # sample a sequence of values for `X_ijk`
    set.seed(seed_val)
    target_x_ijk_samples <- replicate(N_TEST_SAMPLES,
                                      sample_x_ijk(test_data_prior_prob, test_data_posterior_prob))

    # sample a sequence of values for sex on the day before the fertile window
    set.seed(seed_val)
    target_day_before_samples <- as.integer(runif(N_TEST_SAMPLES) < cohort_sex_prob)

    test_data <- c(miss_day_idx   = test_data_miss_day_idx - 1L,
                   day_idx        = test_data_nonmiss_day_idx - 1L,
                   prior_prob     = test_data_prior_prob,
                   posterior_prob = test_data_posterior_prob,
                   xi_i           = test_data_xi_i)

    target_output <- list(x_samples           = X,
                          x_ijk_samples       = target_x_ijk_samples,
                          day_before_samples  = target_day_before_samples,
                          prior_prob_no_prev  = target_prior_prob_no_prev,
                          prior_prob_yes_prev = target_prior_prob_yes_prev,
                          posterior_prob      = target_posterior_prob)

    list(test_data     = test_data,
         target_output = target_output)
}
