get_cyc_idx_list <- function(comb_dat, var_nm) {

    out <- lapply(unique(comb_dat[[var_nm$id]]), function(x) {

        subj_idx <- which(comb_dat[[var_nm$id]] == x)

        lapply(unique(comb_dat[subj_idx, var_nm$cyc]), function(y) {
            subj_idx[comb_dat[subj_idx, var_nm$cyc] == y]
        })
    })

    unlist(out, recursive = FALSE)
}




# takes an atomic vector of type logical, numerical, factor, or character and
# maps it to a logical vector.  Missing values are preserved.

map_vec_to_bool <- function(sex_vec) {

    allowed_yes <- c("y", "Y", "yes", "Yes", "YES")

    if (is.logical(sex_vec)) {
        return(sex_vec)
    }
    else if (is.numeric(sex_vec)) {
        return(sex_vec != 0)
    }
    else if (is.factor(sex_vec) || is.character(sex_vec)) {
        out <- rep(NA, length(sex_vec))
        out[! is.na(sex_vec)] <- sex_vec[! is.na(sex_vec)] %in% allowed_yes
        return(out)
    }
    else {
        # illegal form of intercourse status
        stop("invalid form of sex_vec")
    }
}
