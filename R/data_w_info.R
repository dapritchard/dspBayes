# TODO: this is a half-unfished rewrite of `get_w_day_blocks` (and possible
# merge in `get_w_to_days_idx` and `get_w_cyc_to_subj_idx`?)

# get_w_info <- function() {

#     cyc_idx_list <- get_cyc_idx_list(comb_dat, var_nm)

#     sex_or_miss_bool <- map_vec_to_bool(comb_dat[[var_nm$sex]])
#     is.na(sex_yes_or_miss_bool) <- TRUE

#     preg_bool <- map_vec_to_bool(comb_dat[[var_nm$preg]])

#     # containers for the data we are constructing
#     w_day_blocks <- vector("list", length(cyc_idx_list))
#     w_to_days_idx <- vector("list", length(cyc_idx_list))

#     ctr <- 1L
#     for (curr_cyc_idx in cyc_idx_list) {

#         if (! preg_bool[curr_cyc_idx[1L]]) {
#             next
#         }

#         # else: this is a pregnancy cycle

#         sex_idx <- curr_cyc_idx[sex_or_miss_bool[curr_cyc_idx]]

#         # obtain the needed values for the current block
#         beg_idx <- sex_idx[1L]
#         n_days <- length(curr_idx)
#         subj_idx <- which(unique(comb_dat[[var_nm$id]]) == curr_id)

#         ctr <- ctr + 1L
#     }
# }
