utest_u_categ <- function(dsp_data, u_categ_seed) {

    # TODO
    calc_posterior_u <- function() {
        numeric(0)
    }


    calc_alt_utau_list <- function(utau, tau_coefs, sex_coef, miss_block, cov_miss_x_idx, var_info) {

        block_ptr_idx <- seq(miss_block["beg_sex_idx"] + 1L,
                             miss_block["beg_sex_idx"] + miss_block["n_sex_days"],
                             1L)
        block_cov_miss_x_idx <- cov_miss_x_idx[block_ptr_idx] + 1L
        block_utau <- utau[block_cov_miss_x_idx]

        # prev_sex <- sapply(x_miss_day[block_cov_miss_x_idx], function(z) z["prev"])
        # block_utau_w_sex <- ifelse((prev_sex == 0) | (prev_sex == -2),
        #                            block_utau,
        #                            block_utau + sex_coef)

        tau_star <- ifelse(miss_block["u_col"] == var_info["ref_col"],
                           0,
                           tau_coefs[miss_block["u_col"] + 1L])

        col_idx <- (var_info["col_start"] + 1L) : var_info["ref_col"]
        alt_utau <- lapply(tau_coefs[col_idx], function(tau_j) block_utau + (tau_j - tau_star))
        if (var_info["ref_col"] != var_info["col_end"]) {
            alt_utau$reference <- block_utau - tau_star
        }

        alt_utau
    }


    calc_posterior_x <- function(X, x_miss_day, alt_utau_list, miss_block, cov_miss_x_idx) {

        if (miss_block["beg_sex_idx"] == -1L) {
            return(rep(1, length(alt_utau_list)))
        }

        block_ptr_idx <- seq(miss_block["beg_sex_idx"] + 1L,
                             miss_block["beg_sex_idx"] + miss_block["n_sex_days"],
                             1L)
        block_cov_miss_x_idx <- cov_miss_x_idx[block_ptr_idx] + 1L
        day_idx <- sapply(x_miss_day[block_cov_miss_x_idx], function(z) z["idx"]) + 1L

        prev_sex <- sapply(x_miss_day[block_cov_miss_x_idx], function(z) z["prev"])
        sex_coef_vec <- ifelse((prev_sex == 0) | (prev_sex == -2), 0, sex_coef)

        posterior_x <- vector("numeric", length(alt_utau_list))
        for (j in seq_along(alt_utau_list)) {

            curr_utau <- alt_utau_list[[j]] + sex_coef_vec

            curr_posterior_x_days <- ifelse(X[day_idx] == 1,
                                            1 / (1 + exp(- curr_utau)),
                                            1 / (1 + exp(curr_utau)))
            posterior_x[j] <- prod(curr_posterior_x_days)
        }

        posterior_x
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

    # bind variables to global U data
    X <- dsp_data$intercourse$X
    x_miss_day <- dsp_data$intercourse$miss_day
    u_sex_map <- dsp_data$cov_miss_x_idx
    utau <- dsp_data$utau
    tau_coefs <- dsp_data$tau_fit$u_coefs
    sex_coef <- dsp_data$tau_fit$sex_coef

    # bind variables to variable-specific data
    var_info       <- dsp_data$u_miss_info[[var_idx]]$var_info
    u_prior_probs  <- dsp_data$u_miss_info[[var_idx]]$u_prior_probs
    var_block_list <- dsp_data$u_miss_info[[var_idx]]$var_block_list

    # choose a block with some missing intercourse data, if possible, to use for
    # testing
    beg_sex_idx <- sapply(var_block_list, function(x) x["beg_sex_idx"])
    block_idx <- which(beg_sex_idx >= 0L) %>% head(., 1L) %>% unname
    miss_block <- var_block_list[[ block_idx ]]

    target_data <- c(var_info, c(var_idx   = var_idx - 1L,
                                 n_days    = nrow(dsp_data$U),
                                 block_idx = block_idx - 1L))

    # calculate P(W | U)
    target_posterior_u <- calc_posterior_u()

    # calculate P(X | U)
    alt_utau_list <- calc_alt_utau_list(utau, tau_coefs, sex_coef, miss_block, u_sex_map, var_info)
    target_posterior_x <- calc_posterior_x(X, x_miss_day, alt_utau_list, miss_block, u_sex_map)

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
                           alt_utau_vals = alt_utau_list %>% unlist)

    list(input   = input,
         data    = target_data,
         samples = target_samples)
}
