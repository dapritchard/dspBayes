# just a stub for now
get_gam_info <- function(var_categ_status) {

    var_info <- vector("list", length(var_categ_status))
    for (k in seq_along(var_info)) {
        var_info[[k]] <- list(a=1, b=1, p=0.5, bndL=0, bndU=Inf, is_categ=var_categ_status[k])
    }

    return(var_info)
}


# just a stub for now
get_hyp_phi <- function() {
    list(c1=1, c2=1)
}
