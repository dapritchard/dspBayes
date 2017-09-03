utest_u_categ <- function(dsp_data, u_categ_seed) {

    # TODO
    calc_posterior_u <- function() {
        numeric(0)
    }

    calc_posterior_x <- function(X, tau, alt_utau_vals, miss_block) {
        # stub
        vector("numeric", length(10))
    }

    sample_covariate <- function(posterior_w_probs, posterior_x_probs, u_prior_probs) {

        unnormalized_probs <- vector("numeric", length(u_prior_probs))
        for (i in seq_along(unnormalized_probs)) {
            unnormalized_probs[i] <- (posterior_w_probs[i] *
                                      posterior_x_probs[i] *
                                      u_prior_probs[i])
        }
        probs <- unnormalized_probs / sum(unnormalized_probs)

        u <- runif(1L)

        sum_probs <- probs[1L]
        ctr <- 1L
        while (u > sum_probs) {
            ctr <- ctr + 1L
            sum_probs <- sum_probs + probs[ctr]
        }
        # subtract 1 for 0-based indexing
        ctr - 1L
    }

    # pick a categorical variable to test
    if (! any(dsp_data$u_miss_type == 1L)) {
        stop("need to figure out how to handle no missing categorical vars in testing")
    }
    else {
        var_idx <- which(dsp_data$u_miss_type == 1L) %>% head(., 1L) %>% unname
    }

    # bind variables to variable-specific data
    var_info       <- dsp_data$u_miss_info[[var_idx]]$var_info
    u_prior_probs  <- dsp_data$u_miss_info[[var_idx]]$u_prior_probs
    var_block_list <- dsp_data$u_miss_info[[var_idx]]$var_block_list

    target_data <- c(var_info, c(var_idx   = var_idx - 1L,
                                 n_days    = nrow(dsp_data$U),
                                 block_idx = 0L))                  # TODO: change this idx?

    # calculate P(W | U)
    target_posterior_u <- calc_posterior_u()

    # calculate P(X | U)
    # alt_utau_list <- calc_alt_utau_list()
    target_posterior_x <- calc_posterior_x()

    # a clumsy way to contruct full conditional likelihood values for the W and
    # X.  Note that these need not sum to 1, although they do here.
    z <- seq_along(u_prior_probs)
    input_posterior_w_probs <- z / sum(z)
    input_posterior_x_probs <- rev(z) / sum(z)

    # perform 20 samples for fixed values of P(W | U) and P(X | U) (20 chosen
    # arbitrarily)
    set.seed(u_categ_seed)
    target_sample_covs <- replicate(20L, sample_covariate(input_posterior_w_probs,
                                                          input_posterior_x_probs,
                                                          u_prior_probs))

    input <- list(posterior_w_probs = input_posterior_w_probs,
                  posterior_x_probs = input_posterior_x_probs)

    target_samples <- list(posterior_u   = target_posterior_u,
                           posterior_x   = target_posterior_x,
                           sample_covs   = target_sample_covs,
                           alt_utau_vals = numeric(10))

    list(input   = input,
         data    = target_data,
         samples = target_samples)
}
