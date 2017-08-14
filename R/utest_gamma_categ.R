utest_gamma_categ <- function(dsp_data, W, xi, ubeta, seed_val) {

    log_dgamma_norm_const <- function(a, b) {
        a * log(b) - lgamma(a)
    }

    log_dgamma_trunc_const <- function(a, b, bnd_l, bnd_u) {
        log(pgamma(bnd_u, a, b) - pgamma(bnd_l, a, b))
    }

    calc_log_d2_const_terms <- function(hyp_a, hyp_b, hyp_p, bnd_l, bnd_u) {
        (log(1 - hyp_p) +
         log_dgamma_norm_const(hyp_a, hyp_b) -
         hyp_b -
         log_dgamma_trunc_const(hyp_a, hyp_b, bnd_l, bnd_u))
    }

    calc_a_tilde <- function(dsp_data, W, hyp_a, h) {
        # add 1 because `h` and `w_to_days_idx` use 0-based indexing
        U_h <- dsp_data$U[, h + 1L]
        hyp_a + sum(U_h[dsp_data$w_to_days_idx + 1L] * W)
    }

    calc_b_tilde <- function(dsp_data, hyp_b, xi_day, ubeta_no_h, h) {
        U_h <- dsp_data$U[, h + 1L]
        X <- dsp_data$intercourse$X
        hyp_b + sum(U_h * X * (exp(log(xi_day) + ubeta_no_h)))
    }

    calc_p_tilde <- function(hyp_a, hyp_b, hyp_p, bnd_l, bnd_u, a_tilde, b_tilde) {
        log_d2 <- (log(1 - hyp_p) +
                   log_dgamma_norm_const(hyp_a, hyp_b) +
                   log_dgamma_trunc_const(a_tilde, b_tilde, bnd_l, bnd_u) +
                   (b_tilde - hyp_b) -
                   log_dgamma_norm_const(a_tilde, b_tilde) -
                   log_dgamma_trunc_const(hyp_a, hyp_b, bnd_l, bnd_u))
        hyp_p / (hyp_p + exp(log_d2))
    }

    rgamma_posterior <- function(n, a_tilde, b_tilde, p_tilde, bnd_l, bnd_u) {

        vals <- numeric(n)
        for (i in seq_len(n)) {
            if (runif(1, 0, 1) < p_tilde) {
                vals[i] <- 1
            }
            else {

                if ((bnd_l == 0) && (bnd_u == Inf)) {
                    vals[i] <- rgamma(1, a_tilde, b_tilde)
                }
                else {
                    u <- runif(1, pgamma(bnd_l, a_tilde, b_tilde), pgamma(bnd_u, a_tilde, b_tilde))
                    vals[i] <- qgamma(u, a_tilde, b_tilde)
                }
            }
        }

        vals
    }

    calc_xi_day <- function(xi, subj_day_blocks) {
        n_subj_days <- sapply(subj_day_blocks, function(x) x["n_days"])
        idx <- rep(seq_along(subj_day_blocks), n_subj_days)
        xi[idx]
    }

    # arbitrary choices of data
    beta_prev <- 0.01
    n_test <- 20L
    hyp_a <- 1
    hyp_b <- 2
    hyp_p <- 0.5
    h <- 0

    # expand xi to match daily data
    xi_day <- calc_xi_day(xi, dsp_data$subj_day_blocks)

    # calculate the constant term in the calculation of `p_tilde`
    m_log_d2_const_terms_all <- calc_log_d2_const_terms(hyp_a, hyp_b, hyp_p, 0, Inf)
    m_log_d2_const_terms_zero_one <- calc_log_d2_const_terms(hyp_a, hyp_b, hyp_p, 0, 1)
    m_log_d2_const_terms_one_inf <- calc_log_d2_const_terms(hyp_a, hyp_b, hyp_p, 1, Inf)
    m_log_d2_const_terms_zero_half <- calc_log_d2_const_terms(hyp_a, hyp_b, 0, 0, 0.5)

    # calculate the log of the constant term in a gamma distribution
    log_dgamma_norm_const_val <- log_dgamma_norm_const(hyp_a, hyp_b)

    # calculate the inverse norming constant incurred from truncating a gamma distribution
    log_dgamma_trunc_const_all <- log_dgamma_trunc_const(hyp_a, hyp_b, 0, Inf)
    log_dgamma_trunc_const_zero_one <- log_dgamma_trunc_const(hyp_a, hyp_b, 0, 1)
    log_dgamma_trunc_const_one_inf <- log_dgamma_trunc_const(hyp_a, hyp_b, 1, Inf)
    log_dgamma_trunc_const_zero_half <- log_dgamma_trunc_const(hyp_a, hyp_b, 0, 0.5)

    # calculate `a_tilde`
    a_tilde <- calc_a_tilde(dsp_data, W, hyp_a, h)

    # calculate `b_tilde`
    ubeta_no_h <- ubeta - (dsp_data$U[, h + 1L] * beta_prev)
    b_tilde <- calc_b_tilde(dsp_data, hyp_b, xi_day, ubeta_no_h, h)

    # calculate `p_tilde`
    p_tilde_all <- calc_p_tilde(hyp_a, hyp_b, hyp_p, 0, Inf, a_tilde, b_tilde)
    p_tilde_zero_one <- calc_p_tilde(hyp_a, hyp_b, hyp_p, 0, 1, a_tilde, b_tilde)
    p_tilde_one_inf <- calc_p_tilde(hyp_a, hyp_b, hyp_p, 1, Inf, a_tilde, b_tilde)
    p_tilde_zero_half <- calc_p_tilde(hyp_a, hyp_b, 0, 0, 0.5, a_tilde, b_tilde)

    # calculate a sequence of samples from the posterior distribution
    set.seed(seed_val)
    target_samples_all <- rgamma_posterior(n_test, a_tilde, b_tilde, p_tilde_all, 0, Inf)
    set.seed(seed_val)
    target_samples_zero_one <- rgamma_posterior(n_test, a_tilde, b_tilde, p_tilde_zero_one, 0, 1)
    set.seed(seed_val)
    target_samples_one_inf <- rgamma_posterior(n_test, a_tilde, b_tilde, p_tilde_one_inf, 1, Inf)
    set.seed(seed_val)
    target_samples_zero_half <- rgamma_posterior(n_test, a_tilde, b_tilde, p_tilde_zero_half, 0, 0.5)

    # calculate updated values of `U * beta` based upon first sample in sequence
    ubeta_all       <- ubeta_no_h + (dsp_data$U[, h] * target_samples_all[1L])
    ubeta_zero_one  <- ubeta_no_h + (dsp_data$U[, h] * target_samples_zero_one[1L])
    ubeta_one_inf   <- ubeta_no_h + (dsp_data$U[, h] * target_samples_one_inf[1L])
    ubeta_zero_half <- ubeta_no_h + (dsp_data$U[, h] * target_samples_zero_half[1L])

    # data input for testing routines
    target_data <- c(beta_prev                        = beta_prev,
                     hyp_a                            = hyp_a,
                     hyp_b                            = hyp_b,
                     m_log_d2_const_terms_all         = m_log_d2_const_terms_all,
                     m_log_d2_const_terms_zero_one    = m_log_d2_const_terms_zero_one,
                     m_log_d2_const_terms_one_inf     = m_log_d2_const_terms_one_inf,
                     m_log_d2_const_terms_zero_half   = m_log_d2_const_terms_zero_half,
                     log_dgamma_norm_const_val        = log_dgamma_norm_const_val,
                     log_dgamma_trunc_const_all       = log_dgamma_trunc_const_all,
                     log_dgamma_trunc_const_zero_one  = log_dgamma_trunc_const_zero_one,
                     log_dgamma_trunc_const_one_inf   = log_dgamma_trunc_const_one_inf,
                     log_dgamma_trunc_const_zero_half = log_dgamma_trunc_const_zero_half,
                     a_tilde                          = a_tilde,
                     b_tilde                          = b_tilde,
                     p_tilde_all                      = p_tilde_all,
                     p_tilde_zero_one                 = p_tilde_zero_one,
                     p_tilde_one_inf                  = p_tilde_one_inf,
                     p_tilde_zero_half                = p_tilde_zero_half)

    # gamma specifications
    gamma_specs_all <- c(type  = 0,
                         h     = h,
                         hyp_a = hyp_a,
                         hyp_b = hyp_b,
                         hyp_p = hyp_p,
                         bnd_l = 0,
                         bnd_u = Inf)
    gamma_specs_zero_one <- replace(gamma_specs_all, "bnd_u", 1)
    gamma_specs_one_inf <- replace(gamma_specs_all, "bnd_l", 1)
    gamma_specs_zero_half <- replace(gamma_specs_all, c("hyp_p", "bnd_u"), c(0, 0.5))
    input_specs <- list(gamma_specs_all       = gamma_specs_all,
                        gamma_specs_zero_one  = gamma_specs_zero_one,
                        gamma_specs_one_inf   = gamma_specs_one_inf,
                        gamma_specs_zero_half = gamma_specs_zero_half)

    # target output for sample generation
    target_samples <- list(target_samples_all       = target_samples_all,
                           target_samples_zero_one  = target_samples_zero_one,
                           target_samples_one_inf   = target_samples_one_inf,
                           target_samples_zero_half = target_samples_zero_half,
                           target_ubeta_no_h        = ubeta_no_h,
                           target_ubeta_all         = ubeta_all,
                           target_ubeta_zero_one    = ubeta_zero_one,
                           target_ubeta_one_inf     = ubeta_one_inf,
                           target_ubeta_zero_half   = ubeta_zero_half)

    list(input_specs    = input_specs,
         target_data    = target_data,
         target_samples = target_samples)
}
