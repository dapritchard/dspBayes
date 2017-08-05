extract_var_nm <- function(dsp_model,
                           baseline,
                           cycle,
                           daily,
                           id_name,
                           cyc_name,
                           sex_name,
                           fw_name) {

    # `all_nm` is the union of all of the variables that we need across the
    # various datasets.  `all_base`, `all_cyc`, and `all_day` are the names of
    # all of the variables in the model that are in the various datasets.  Note
    # that these expressions all handle NULL inputs without choking.
    all_nm <- c(id_name, cyc_name, sex_name, fw_name, all.vars(dsp_model)) %>% unique
    all_base <- colnames(baseline)[colnames(baseline) %in% all_nm]
    all_cyc <- colnames(cycle)[colnames(cycle) %in% all_nm]
    all_day <- colnames(daily)[colnames(daily) %in% all_nm]

    # these can be NULL if the object or the column names are NULL.  Convert to
    # an empty vector for consistency.
    if (is.null(all_base)) all_base <- character(0L)
    if (is.null(all_cyc)) all_cyc <- character(0L)
    if (is.null(all_day)) all_day <- character(0L)

    # collect into a single list and ensure that the variable names make sense
    var_nm <- list(id       = id_name,
                   cyc      = cyc_name,
                   sex      = sex_name,
                   fw       = fw_name,
                   preg     = all.vars(dsp_model)[1L],
                   all      = all_nm,
                   all_base = all_base,
                   all_cyc  = all_cyc,
                   all_day  = all_day)
    check_var_nm(baseline, cycle, daily, var_nm)

    var_nm
}




check_var_nm <- function(baseline, cycle, daily, var_nm) {

    # check that join variables exist for all datasets
    if (! is.null(baseline) && !(var_nm$id %in% colnames(baseline))) {
        stop("ID variable ", var_nm$id, " not found in baseline data", call. = FALSE)
    }
    if (! is.null(cycle) && !(var_nm$id %in% colnames(cycle))) {
        stop("ID variable ", var_nm$id, " not found in cycle data", call. = FALSE)
    }
    if (! is.null(cycle) && !(var_nm$cyc %in% colnames(cycle))) {
        stop("cycle variable ", var_nm$cyc, " not found in cycle data", call. = FALSE)
    }
    if (! (var_nm$id %in% colnames(daily))) {
        stop("ID variable ", var_nm$id, " not found in daily data", call. = FALSE)
    }
    if (! (var_nm$cyc %in% colnames(daily))) {
        stop("ID variable ", var_nm$cyc, " not found in daily data", call. = FALSE)
    }

    # daily data must have both intercourse data and fertile window data
    if (! (var_nm$sex %in% colnames(daily))) {
        stop("intercourse variable ", var_nm$sex, " not found in daily data", call. = FALSE)
    }
    if (! (var_nm$fw %in% colnames(daily))) {
        stop("fertile window variable ", var_nm$id, " not found in daily data", call. = FALSE)
    }

    # pregnancy variable must either be in cycle or daily
    if ((! is.null(cycle) || !(var_nm$preg %in% cycle)) && !(var_nm$preg %in% daily)) {
        stop("pregnancy variable ", var_nm$preg,
             " not found in either in cycle data daily data or in daily data", call. = FALSE)
    }

    # variables other than the join variables for the various datasets.  NULL
    # detasets will have a character(0L) returned for this.
    non_join_base <- setdiff(var_nm$all_base, var_nm$id)
    non_join_cyc <- setdiff(var_nm$all_cyc, c(var_nm$id, var_nm$cyc))
    non_join_day <- setdiff(var_nm$all_day, c(var_nm$id, var_nm$cyc))

    # it is ambiguous to have the same (non-join) variable in the model in
    # multiple datasets
    if (! identical(intersect(non_join_base, non_join_cyc), character(0L))) {
        stop("baseline and cycle data have the same variable(s) in the model: ",
             paste(intersect(non_join_base, non_join_cyc), collapse = ", "))
    }
    if (! identical(intersect(non_join_base, non_join_day), character(0L))) {
        stop("baseline and daily data have the same variable(s) in the model: ",
             paste(intersect(non_join_base, non_join_day), collapse = ", "))
    }
    if (! identical(intersect(non_join_cyc, non_join_day), character(0L))) {
        stop("cycle and daily data have the same variable(s) in the model: ",
             paste(intersect(non_join_cyc, non_join_day), collapse = ", "))
    }

    # if there are variables in the model that are not in the datasets then this
    # will cause a crash later
    union_vars <- c(var_nm$base, var_nm$cyc, var_nm$day) %>% unique
    if (! isTRUE(all.equal(union_vars, var_nm$all, check.attributes = FALSE))) {
        stop("There are variables in the model that are not in any of the datasets:  ",
             paste(setdiff(var_nm$all, union_vars), collapse = ", "))
    }
}
