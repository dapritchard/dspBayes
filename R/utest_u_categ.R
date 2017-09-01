utest_u_categ <- function(dsp_data) {

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


    target_data_u <- c(var_info, c(var_idx = var_idx,
                                   n_days  = nrow(dsp_data$U)))


    list(target_data_u = target_data_u)
}
