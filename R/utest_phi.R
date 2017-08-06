utest_phi <- function(xi_vals,
                      hyp_c1       = 1,
                      hyp_c2       = 1,
                      phi_proposal = 1.05,
                      phi_init     = hyp_c1 / hyp_c2,
                      delta        = 0.1,
                      n_phi_samp   = 20L,
                      seed_val     = 99L) {

    # sample phi from proposal dist
    samp_phi_prop <- function(phi, delta, prop_dist_nm = "uniform") {
        if (prop_dist_nm == "normal") {
            rnorm(n = 1, mean = phi, sd = delta) %>% abs
        } else {
            runif(n = 1, min = phi - delta, max = phi + delta) %>% abs
        }
    }

    # calculate log(r), where r is given by:
    #
    #     ( prod_i pi(xi_i | proposal_val) ) * pi(proposal_val)
    #     -----------------------------------------------------
    #         ( prod_i pi(xi_i | curr_val) ) * pi(curr_val)

    calc_log_r <- function(xi, curr_phi, proposal_val, hyp_c1, hyp_c2) {

        logNumerTerm1 <- dgamma(xi, shape = proposal_val, rate = proposal_val, log = TRUE) %>% sum
        logNumerTerm2 <- dgamma(proposal_val, shape = hyp_c1, rate = hyp_c2, log = TRUE)

        logDenomTerm1 <- dgamma(xi, shape = curr_phi, rate = curr_phi, log = TRUE ) %>% sum
        logDenomTerm2 <- dgamma(curr_phi, shape = hyp_c1, rate = hyp_c2, log = TRUE)

        logNumerTerm1 + logNumerTerm2 - logDenomTerm1 - logDenomTerm2
    }

    # sample a new value of phi
    samp_phi <- function(curr_phi, xi, hyp_c1, hyp_c2, delta, prop_dist_nm) {

        proposal_val <- samp_phi_prop(curr_phi, delta)
        log_r <- calc_log_r(xi, curr_phi, proposal_val, hyp_c1, hyp_c2)
        if (log(runif(1L)) < log_r) {
            proposal_val
        } else {
            curr_phi
        }
    }

    # log_dgamma_norm_const
    out_log_dgamma_norm_const <- log(phi_proposal^phi_proposal / gamma(phi_proposal))

    # calc_log_proportion_dgamma_phi
    out_log_proportion_dgamma_phi <-
        (dgamma(phi_proposal, shape = hyp_c1, rate = hyp_c2, log = TRUE) -
         dgamma(phi_init, shape = hyp_c1, rate = hyp_c2, log = TRUE))

    # calc_log_proportion_dgamma_xi
    out_log_proportion_dgamma_xi <-
        (sum( dgamma(xi_vals, shape = phi_proposal, rate = phi_proposal, log = TRUE) ) -
         sum( dgamma(xi_vals, shape = phi_init, rate = phi_init, log = TRUE) ))

    # m_log_norm_const
    out_m_log_norm_const <- log(phi_init^phi_init / gamma(phi_init))

    # calc_log_r
    # hyp_phi_list <- list(c1 = hyp_c1, c2 = hyp_c2)
    out_calc_log_r <- calc_log_r(xi_vals, phi_init, phi_proposal, hyp_c1, hyp_c2)

    # sample a series of phi updates.  Note: this is a random sample, so we have
    # to set a seed.
    set.seed(seed_val)
    target_samples <- vector("numeric", n_phi_samp)
    curr_phi <- phi_init
    accept_ctr <- 0L
    for (i in seq_len(n_phi_samp)) {

        target_samples[i] <- samp_phi(curr_phi, xi_vals, hyp_c1, hyp_c2, delta)

        if (target_samples[i] != curr_phi) {
            accept_ctr <- accept_ctr + 1L
        }
        curr_phi <- target_samples[i]
    }

    # collect data
    target_data <- c(phi_val                   = phi_init,
                     proposal_val              = phi_proposal,
                     log_dgamma_norm_const     = out_log_dgamma_norm_const,
                     log_proportion_dgamma_phi = out_log_proportion_dgamma_phi,
                     log_proportion_dgamma_xi  = out_log_proportion_dgamma_xi,
                     m_log_norm_const          = out_m_log_norm_const,
                     calc_log_r                = out_calc_log_r,
                     accept_ctr                = accept_ctr,
                     seed_val                  = seed_val)

    # return testing inputs and results
    list(target_data    = target_data,
         target_samples = target_samples)
}
