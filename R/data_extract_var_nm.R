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

    # collect into a single list
    list(id       = id_name,
         cyc      = cyc_name,
         sex      = sex_name,
         fw       = fw_name,
         preg     = all.vars(dsp_model)[1L],
         all      = all_nm,
         all_base = all_base,
         all_cyc  = all_cyc,
         all_day  = all_day)
}



# TODO: check that these make sense.
#
# Do each of the non-NULL datasets have variables of interest in them and the
# appropriate join variables?
#
# is there a column for at minimum id, cyc, fw, preg, sex?
#
# maybe allow and just check for conflicts

#Call checks from the above fcn
