get_tau_fit <- function(comb_dat, var_nm, dsp_model, use_na) {

    # case: we're not imputing missing for sex, so return some empty data
    if (! ((use_na == "sex") || (use_na == "all"))) {
        return(list(u_coefs  = numeric(0L),
                    sex_coef = 0))
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

    list(u_coefs  = tau_coefs[-1L],
         sex_coef = tau_coefs[1L])
}




get_utau <- function(U, tau_fit, intercourse_data) {

    # case: we're not imputing missing for sex, so return some empty data
    if (! ((use_na == "sex") || (use_na == "all"))) {
        return(numeric(0L))
    }

    # assert that the tau model fit is consistent with the design matrix
    # expansion of the data, U.  This is an internal check since these should
    # always be consistent unless there is a bug in the code.
    if (! identical(colnames(U), names(tau))) {
        stop('internal inconsistency, the names for "U" do not match those for "tau"\n',
             'U:  ', colnames(U), '\n', 'tau:  ', names(tau))
    }

    # calculates `(U, sex_yester) %*% (tau, sex_coef)`.  The reason that
    # `sex_yester * sex_coef` term is calculated using the
    U %*% tau$u_coefs + ifelse(intercourse_data$sex_yester == 1L, tau$sex_coef, 0)
}
