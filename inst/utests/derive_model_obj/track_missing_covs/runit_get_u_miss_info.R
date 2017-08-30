# construct data ---------------------------------------------------------------

# idx:    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
# ---------------------------------------------------------------------------------------------
id   <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6)
cyc  <- c(1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2)
preg <- c(0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)
sex  <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
age  <- c(2, 2, 2, 2, 2, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 2, 2, 9, 9, 9, 9, 5, 5, 5, 5, 5, 5, 5, 5)
bmi  <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2)
edu  <- c(2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
opk  <- c(0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0)
drnk <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)

# all vars are factors except for `age`
bmi  <- as.factor(bmi)
edu  <- as.factor(edu)
opk  <- as.factor(opk)
drnk <- as.factor(drnk)

# all vars have missing except for `edu`
base_blocks <- list(1:5,      6:8, 9:15,               16:17, 18:21,        22:29)
cyc_blocks <-  list(1:2, 3:5, 6:8, 9:10, 11:13, 14:15, 16:17, 18:19, 20:21, 22:27, 28:29)
sex[c(7, 13, 14, 17, 20, 21, 25)] <- NA_real_
age[base_blocks[1:3] %>% unlist] <- NA_real_
bmi[base_blocks[5:6] %>% unlist] <- NA_integer_
opk[cyc_blocks[c(1, 2, 3, 6, 9)] %>% unlist] <- NA_integer_
drnk[c(9, 12, 14, 18, 20, 26, 27, 29)] <- NA_integer_

comb_dat <- data.frame(id   = id,
                       cyc  = cyc,
                       preg = preg,
                       sex  = sex,
                       age  = age,
                       bmi  = bmi,
                       edu  = edu,
                       opk  = opk,
                       drnk = drnk)

model_with_intercept <- formula(~ age + bmi + edu + opk + drnk)
model_no_intercept <- formula(~ 0 + age + bmi + edu + opk + drnk)

expan_with_intercept <- expand_model_rhs(comb_dat, model_with_intercept)
expan_no_intercept <- expand_model_rhs(comb_dat, model_no_intercept)

cov_with_intercept_info <- get_cov_col_miss_info(expan_with_intercept, model_with_intercept, "all")
cov_no_intercept_info <- get_cov_col_miss_info(expan_no_intercept, model_no_intercept, "all")

var_nm <- list(id       = "id",
               cyc      = "cyc",
               preg     = "preg",
               sex      = "sex",
               all_base = c("id", "age", "bmi"),
               all_cyc  = c("id", "cyc", "preg", "opk"))

# TODO: remove data duplication above (below is new)

row_info_with_intercept <- get_cov_row_miss_info(comb_dat,
                                                 var_nm,
                                                 expan_with_intercept,
                                                 cov_with_intercept_info)

row_info_no_intercept <- get_cov_row_miss_info(comb_dat,
                                               var_nm,
                                               expan_no_intercept,
                                               cov_no_intercept_info)




# target output ----------------------------------------------------------------

var_info_nm <- c("col_start", "col_end", "ref_col", "n_categs", "max_n_days_miss", "max_n_sex_days_miss")

max_miss <- sapply(row_info_with_intercept$cov_row_miss_list,
                   function(x) sapply(x, function(y) y["n_days"]) %>% max)

var_info_with_intercept <- list(age  = setNames(c(1L, 2L, 2L, 1L, 7L, 2L), var_info_nm),
                                bmi  = setNames(c(2L, 6L, 5L, 4L, 8L, 2L), var_info_nm),
                                opk  = setNames(c(7L, 9L, 8L, 2L, 3L, 2L), var_info_nm),
                                drnk = setNames(c(8L,10L, 9L, 2L, 1L, 1L), var_info_nm))

var_info_no_intercept <- list(age  = setNames(c(0L, 1L, 1L, 1L, 7L, 2L), var_info_nm),
                              bmi  = setNames(c(1L, 5L, 5L, 4L, 8L, 2L), var_info_nm),
                              opk  = setNames(c(7L, 9L, 8L, 2L, 3L, 2L), var_info_nm),
                              drnk = setNames(c(8L,10L, 9L, 2L, 1L, 1L), var_info_nm))

target_with_intercept <- target_no_intercept <- vector("list", 4L)

for (i in seq_along(target_with_intercept)) {

    target_with_intercept[[i]] <- list(var_info       = var_info_with_intercept[[i]],
                                       u_prior_probs  = cov_with_intercept_info[[i]]$empirical_probs,
                                       var_block_list = row_info_with_intercept$cov_row_miss_list[[i]])

    target_no_intercept[[i]] <- list(var_info       = var_info_no_intercept[[i]],
                                     u_prior_probs  = cov_no_intercept_info[[i]]$empirical_probs,
                                     var_block_list = row_info_no_intercept$cov_row_miss_list[[i]])
}




# begin testing ----------------------------------------------------------------

test_get_u_miss_info <- function() {

    out_with_intercept <- get_u_miss_info(cov_with_intercept_info, row_info_with_intercept)
    identical(target_with_intercept, out_with_intercept)

    out_no_intercept <- get_u_miss_info(cov_no_intercept_info, row_info_no_intercept)
    identical(target_no_intercept, out_no_intercept)
}
