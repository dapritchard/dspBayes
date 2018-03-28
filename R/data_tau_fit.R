get_tau_fit <- function(comb_dat, var_nm, dsp_model, use_na) {

    # case: we're not imputing missing for sex, so return some empty data
    if (! ((use_na == "sex") || (use_na == "all"))) {
        return(list(u_coefs         = numeric(0L),
                    sex_coef        = 0,
                    cohort_sex_prob = 0))
    }

    # strips the response variable (i.e. pregnancy) from the model, replaces it
    # with the intercourse variable, and adds sex yesterday status as a
    # predictor
    stripped_response <- delete.response(terms(dsp_model))
    sex_outcome <- reformulate(termlabels = attr(stripped_response, "term.labels"),
                               response   = var_nm$sex,
                               intercept  = attr(stripped_response, "intercept"))
    sex_model <- update(sex_outcome, ~ sex_yester + .)

    # manually transform sex and sex yesterday variables to ensure consistency
    # with our interpretation of the data
    comb_dat[[var_nm$sex]] <- map_vec_to_bool(comb_dat[[var_nm$sex]]) %>% as.integer
    comb_dat[["sex_yester"]] <- map_vec_to_bool(comb_dat[["sex_yester"]]) %>% as.integer

    # fit the model and extract the variable coefficients
    sex_outcome_fit <- glm(sex_model, family = "binomial", data = comb_dat)
    tau_coefs <- coef(sex_outcome_fit)
    # TODO: check if fit is the same as for geeglm? should be using this?
    # TODO: we have to think about what to do if most or all vars have at least 1 missing

    # overall proportion of intercourse in the cohort
    cohort_sex_prob <- mean(comb_dat[[var_nm$sex]], na.rm = TRUE)

    list(u_coefs         = tau_coefs[-1L],
         sex_coef        = tau_coefs[1L],
         cohort_sex_prob = cohort_sex_prob)
}




get_utau <- function(U, tau_fit, xmiss, use_na) {

    # case: we're not imputing missing for sex, so return some empty data
    if (((use_na != "sex") && (use_na != "all")) || (length(xmiss) == 0L)) {
        return(numeric(0L))
    }

    # assert that the tau model fit is consistent with the design matrix
    # expansion of the data, U.  This is an internal check since these should
    # always be consistent unless there is a bug in the code.
    if (! identical(colnames(U), names(tau_fit$u_coefs))) {
        stop('internal inconsistency, the names for "U" do not match those for "tau_fit"\n',
             'U:  ', colnames(U), '\ntau_fit:  ', names(tau_fit$u_coefs))
    }

    # calculate an initial `U * tau` vector for the rows with a missing value of
    # `X`.  The conditional `xmiss >= 1` has the effect of testing whether there
    # is a missing value, since these values are indices in the data, while
    # yes/no sex are denoted by nonpositive values.
    tau <- matrix(tau_fit$u_coefs, ncol = 1L)
    miss_idx <- xmiss[xmiss >= 1L]
    drop(U[miss_idx, , drop = FALSE] %*% tau)
}




# get_utau_old <- function(U, tau_fit, intercourse_data, use_na) {

#     # case: we're not imputing missing for sex, so return some empty data
#     if (((use_na != "sex") && (use_na != "all")) || (length(intercourse_data$miss_day) == 0L)) {
#         return(numeric(0L))
#     }

#     # assert that the tau model fit is consistent with the design matrix
#     # expansion of the data, U.  This is an internal check since these should
#     # always be consistent unless there is a bug in the code.
#     if (! identical(colnames(U), names(tau_fit$u_coefs))) {
#         stop('internal inconsistency, the names for "U" do not match those for "tau_fit"\n',
#              'U:  ', colnames(U), '\ntau_fit:  ', names(tau_fit$u_coefs))
#     }

#     tau <- matrix(tau_fit$u_coefs, ncol = 1L)
#     miss_idx <- sapply(intercourse_data$miss_day, function(x) x["idx"] + 1L)
#     drop(U[miss_idx, , drop = FALSE] %*% tau)
# }
