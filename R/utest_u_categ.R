utest_u_categ <- function(dsp_data, W, xi, beta_coefs, X, u_categ_seed) {

    sample <- function(U,
                       W,
                       xi,
                       var_info,
                       u_prior_probs,
                       var_block_list,
                       ubeta,
                       beta_coefs,
                       utau,
                       tau_coefs,
                       sex_coef,
                       cov_miss_w_idx,
                       cov_miss_x_idx) {

        categ_update <- vector("integer", length(var_block_list))
        ctr <- 1L

        for (miss_block in var_block_list) {

            alt_exp_ubeta_list <- calc_alt_exp_ubeta_list(ubeta, beta_coefs, miss_block, var_info)
            posterior_w <- calc_posterior_w(W, xi, X, alt_exp_ubeta_list, miss_block, cov_miss_w_idx)

            if (miss_block["n_sex_days"] == 0L) {
                posterior_x <- rep(1, length(posterior_w))
            }
            else{
                alt_utau_list <- calc_alt_utau_list(utau,
                                                    tau_coefs,
                                                    sex_coef,
                                                    miss_block,
                                                    cov_miss_x_idx,
                                                    var_info)
                posterior_x <- calc_posterior_x(X,
                                                x_miss_day,
                                                alt_utau_list,
                                                miss_block,
                                                cov_miss_x_idx)
            }

            u_categ <- sample_covariate(posterior_w, posterior_x, u_prior_probs)
            u_col <- u_categ + var_info["col_start"]

            if (u_col != miss_block["u_col"]) {

                # update U
                row_idx <- seq(miss_block["beg_day_idx"] + 1L,
                               miss_block["beg_day_idx"] + miss_block["n_days"],
                               1)
                # zero out old column
                if (miss_block["u_col"] != var_info["ref_col"]) {
                    U[row_idx, miss_block["u_col"] + 1L] <- 0
                }
                # place 1's in new column
                if (u_col != var_info["ref_col"]) {
                    U[row_idx, u_col + 1L] <- 1
                }

                # update ubeta
                ubeta[row_idx] <- alt_exp_ubeta_list[[u_categ + 1L]] %>% log

                # update utau
                if (miss_block["n_sex_days"] > 0L) {
                    block_ptr_idx <- seq(miss_block["beg_sex_idx"] + 1L,
                                         miss_block["beg_sex_idx"] + miss_block["n_sex_days"],
                                         1L)
                    block_cov_miss_x_idx <- cov_miss_x_idx[block_ptr_idx] + 1L
                    utau[block_cov_miss_x_idx] <- alt_utau_list[[u_categ + 1L]]
                }

                # update block
                categ_update[ctr] <- miss_block["u_col"] <- u_col
            }
            else {
                categ_update[ctr] <- miss_block["u_col"]
            }

            ctr <- ctr + 1L
        }

        list(categ_update = categ_update,
             ubeta_update = ubeta,
             utau_update  = utau,
             u_update     = U[, (var_info["col_start"] + 1L) : var_info["ref_col"], drop = FALSE])
    }


    calc_alt_exp_ubeta_list <- function(ubeta, beta_coefs, miss_block, var_info) {

        block_day_idx <- seq(miss_block["beg_day_idx"] + 1L,
                             miss_block["beg_day_idx"] + miss_block["n_days"],
                             1L)
        block_ubeta <- ubeta[block_day_idx]

        beta_star <- ifelse(miss_block["u_col"] == var_info["ref_col"],
                            0,
                            beta_coefs[miss_block["u_col"] + 1L])

        col_idx <- (var_info["col_start"] + 1L) : var_info["ref_col"]
        alt_ubeta <- lapply(beta_coefs[col_idx], function(beta_j) block_ubeta + (beta_j - beta_star))
        if (var_info["ref_col"] != var_info["col_end"]) {
            alt_ubeta$reference <- block_ubeta - beta_star
        }

        lapply(alt_ubeta, exp)
    }


    calc_posterior_w <- function(W, xi, X, alt_exp_ubeta_list, miss_block, cov_miss_w_idx) {

        xi_i <- xi[miss_block["subj_idx"] + 1L]

        # rows in the `X` data
        block_x_idx <- seq(miss_block["beg_day_idx"] + 1L,
                           miss_block["beg_day_idx"] + miss_block["n_days"],
                           1L)
        x_vals <- X[block_x_idx]
        if (all(x_vals == 0L)) {
            return(1)
        }

        # rows in the map to the `W` data
        block_ptr_idx <- seq(miss_block["beg_w_idx"] + 1L,
                             miss_block["beg_w_idx"] + miss_block["n_days"],
                             1L)
        # rows in the W data (or 0 if occurred during a non-pregnancy cycle)
        block_cov_w_idx <- cov_miss_w_idx[block_ptr_idx] + 1L
        preg_cycle_day_bool <- (block_cov_w_idx != 0)

        # each element in the list is posterior value of W with respect to the
        # `j`-th vector of `exp(U * beta)` values
        posterior_w <- vector("numeric", length(alt_exp_ubeta_list))
        for (j in seq_along(alt_exp_ubeta_list)) {

            # if the day occurred during a non-pregnancy cycle then the value of
            # W is 0, otherwise we have to look it up
            w_vals <- vector("numeric", length(block_cov_w_idx))
            w_vals[preg_cycle_day_bool] <- W[block_cov_w_idx[preg_cycle_day_bool]]

            # remove days in which intercourse didn't occur
            w_vals <- w_vals[x_vals == 1L]
            mean_vals <- xi_i * alt_exp_ubeta_list[[j]][x_vals == 1L]

            # calculate log-likelihoods and then sum and exponentiate
            posterior_w[j] <- dpois(w_vals, mean_vals, TRUE) %>% sum %>% exp
        }

        posterior_w
    }


    calc_alt_utau_list <- function(utau, tau_coefs, sex_coef, miss_block, cov_miss_x_idx, var_info) {

        if (miss_block["n_sex_days"] == 0L) {
            return(lapply(seq_len(var_info["n_categs"]), function(x) numeric(0L)))
        }

        block_ptr_idx <- seq(miss_block["beg_sex_idx"] + 1L,
                             miss_block["beg_sex_idx"] + miss_block["n_sex_days"],
                             1L)
        block_cov_miss_x_idx <- cov_miss_x_idx[block_ptr_idx] + 1L
        block_utau <- utau[block_cov_miss_x_idx]

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

        # rows in the "missing U and missing intercourse" map to the missing
        # intercourse data
        block_ptr_idx <- seq(miss_block["beg_sex_idx"] + 1L,
                             miss_block["beg_sex_idx"] + miss_block["n_sex_days"],
                             1L)
        # rows in the missing intercourse data
        block_cov_miss_x_idx <- cov_miss_x_idx[block_ptr_idx] + 1L
        # rows in the full data
        day_idx <- sapply(x_miss_day[block_cov_miss_x_idx], function(z) z["idx"]) + 1L

        # a vector with an element for each day in the missing block with the
        # value of 0 if no sex occurred yesterday or `sex_coef` if yes
        prev_sex <- sapply(x_miss_day[block_cov_miss_x_idx], function(z) z["prev"])
        sex_coef_vec <- ifelse((prev_sex == 0) | (prev_sex == -2), 0, sex_coef)

        # each element in the list is posterior value of X with respect to the
        # `j`-th vector of `U * tau` values
        posterior_x <- vector("numeric", length(alt_utau_list))
        for (j in seq_along(alt_utau_list)) {

            curr_utau <- alt_utau_list[[j]] + sex_coef_vec

            # calculate `P(X = x | U * tau)` for either `x = 1` or `x = 0` as
            # appropriate
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
    if (any(dsp_data$u_miss_type == 1L)) {
        var_idx <- which(dsp_data$u_miss_type == 1L) %>% head(., 1L) %>% unname
    }
    # if there's no missing U then exit early
    else {
        return(list(input   = list(),
                    data    = integer(0L),
                    samples = list()))
    }

    # bind variables to global U data
    X <- dsp_data$intercourse$X
    x_miss_day <- dsp_data$intercourse$miss_day
    u_preg_map <- dsp_data$cov_miss_w_idx
    u_sex_map <- dsp_data$cov_miss_x_idx
    utau <- dsp_data$utau
    tau_coefs <- dsp_data$tau_fit$u_coefs
    sex_coef <- dsp_data$tau_fit$sex_coef

    # bind variables to variable-specific data
    var_info       <- dsp_data$u_miss_info[[var_idx]]$var_info
    u_prior_probs  <- dsp_data$u_miss_info[[var_idx]]$log_u_prior_probs %>% exp
    var_block_list <- dsp_data$u_miss_info[[var_idx]]$var_block_list

    # choose a block with some missing intercourse data, if possible, to use for
    # testing
    # TODO: can we look for one with a pregnancy cycle also?
    beg_sex_idx <- sapply(var_block_list, function(x) x["beg_sex_idx"])
    block_idx <- which(beg_sex_idx >= 0L) %>% head(., 1L) %>% unname
    if (length(block_idx) == 0L) {
        block_idx <- 1L
    }
    miss_block <- var_block_list[[ block_idx ]]

    target_data <- c(var_info, c(var_idx   = var_idx - 1L,
                                 n_days    = nrow(dsp_data$U),
                                 block_idx = block_idx - 1L))

    # calculate P(W | U)
    ubeta <- dsp_data$U %*% as.matrix(beta_coefs) %>% drop
    alt_exp_ubeta_list <- calc_alt_exp_ubeta_list(ubeta, beta_coefs, miss_block, var_info)
    target_posterior_w <- calc_posterior_w(W, xi, X, alt_exp_ubeta_list, miss_block, u_preg_map) %>% log

    # calculate P(X | U)
    alt_utau_list <- calc_alt_utau_list(utau, tau_coefs, sex_coef, miss_block, u_sex_map, var_info)
    target_posterior_x <- calc_posterior_x(X, x_miss_day, alt_utau_list, miss_block, u_sex_map) %>% log

    # a clumsy way to contruct full conditional likelihood values for the W and
    # X.  Note that these need not sum to 1, although they do here.
    z <- seq_along(u_prior_probs)
    input_posterior_w_probs <- z / sum(z)
    input_posterior_x_probs <- rev(z) / sum(z)

    # perform 20 samples for fixed values of P(W | U) and P(X | U) (the value 20
    # is chosen arbitrarily)
    set.seed(u_categ_seed)
    target_sample_covs <- replicate(20L, sample_covariate(input_posterior_w_probs,
                                                          input_posterior_x_probs,
                                                          u_prior_probs))

    # sample covariates
    set.seed(u_categ_seed)
    out_sample <- sample(dsp_data$U,
                         W,
                         xi,
                         var_info,
                         u_prior_probs,
                         var_block_list,
                         ubeta,
                         beta_coefs,
                         utau,
                         tau_coefs,
                         sex_coef,
                         u_preg_map,
                         u_sex_map)

    input <- list(posterior_w_probs = log(input_posterior_w_probs),
                  posterior_x_probs = log(input_posterior_x_probs))

    target_samples <- list(posterior_w        = target_posterior_w,
                           posterior_x        = target_posterior_x,
                           sample_covs        = target_sample_covs,
                           alt_exp_ubeta_vals = alt_exp_ubeta_list %>% unlist,
                           alt_utau_vals      = alt_utau_list %>% unlist,
                           categ_update       = out_sample$categ_update,
                           ubeta_update       = out_sample$ubeta_update,
                           utau_update        = out_sample$utau_update,
                           u_update           = out_sample$u_update)

    list(input   = input,
         data    = target_data,
         samples = target_samples)
}
