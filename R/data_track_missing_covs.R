get_missing_var_info <- function(expanded_df, dsp_model) {

    # the assign attribute of `model.matrix` is an integer vector with an entry
    # for each column in the design matrix, and corresponding to the term in the
    # formula (stored in `explan_var_nm`) which gave rise to the column.  A
    # value of 0 corresponds to the intercept, if any.
    map_expand_to_orig <- attr(expanded_df, "assign")
    unique_map_vals <- setdiff(map_expand_to_orig, 0L)

    # `terms()` provides a representation of the symbolic model that we can use
    # to access some of the information about the model
    dsp_terms <- terms(dsp_model)
    # variable names for the model RHS
    explan_var_nm <-  attr(dsp_terms, "term.labels")

    # categorical data variable names
    categ_nm <- attr(expanded_df, "contrasts") %>% names

    # container in which to place missing covariate information.  Add 1 in
    # case we need to add an intercept later.
    cov_miss_list <- vector("list", length(unique_map_vals))
    names(cov_miss_list) <- explan_var_nm

    # each iteration stores information for the k-th variable in a list
    # providing the indices in the design matrix that the variable corresponds
    # to, whether the variable is categorical, and if the variable has a
    # reference cell coding
    for (k in unique_map_vals) {

        cov_miss_list[[k]] <- list()

        # current variable name and corresponding column indices in the design
        # matrix
        curr_var_nm <- explan_var_nm[k]
        curr_idx <- which(map_expand_to_orig == k)
        curr_obs <- expanded_df[, curr_idx, drop = FALSE]

        # store current variable column indices using 0-based indexing
        cov_miss_list[[k]]$col_start <- head(curr_idx, 1L) - 1L
        cov_miss_list[[k]]$col_end <- tail(curr_idx, 1L) - 1L

        # case: not a categorical variable
        if (! (curr_var_nm %in% categ_nm)) {
            cov_miss_list[[k]]$categ <- FALSE
            cov_miss_list[[k]]$all_cells <- TRUE
        }
        # case: a categorical variable.  Find out if it is also a reference cell
        # coding or not.
        else {

            cov_miss_list[[k]]$categ <- TRUE

            if (all(rowSums(curr_obs) == 1L)) {
                cov_miss_list[[k]]$all_cells <- TRUE
            } else {
                cov_miss_list[[k]]$all_cells <- FALSE
            }
        }

        # store a vector of missing observation indices
        cov_miss_list[[k]]$miss <- which(! complete.cases(curr_obs))
    }

    # add in an element for the intercept if the model includes one
    if (0L %in% map_expand_to_orig) {
        cov_miss_list$intercept <- list()
        cov_miss_list$intercept$col_start <- 0L
        cov_miss_list$intercept$col_end <- 0L
        cov_miss_list$intercept$categ <- TRUE
        cov_miss_list$intercept$all_cells <- TRUE
        cov_miss_list$intercept$miss <- integer(0L)
    }

    cov_miss_list
}




# is_daily <- function(nm, var_nm) {
#     nm %in% var_nm$all_day
# }




# is_cycle <- function(nm, var_nm) {
#     nm %in% var_nm$all_cyc
# }




# is_baseline <- function(nm, var_nm) {
#     nm %in% var_nm$all_base
# }
