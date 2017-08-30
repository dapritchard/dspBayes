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

# TODO: remove data duplication above (below is new)

cov_with_intercept_info <- get_cov_col_miss_info(expan_with_intercept, model_with_intercept, "all")
cov_no_intercept_info <- get_cov_col_miss_info(expan_no_intercept, model_no_intercept, "all")

var_nm <- list(id       = "id",
               cyc      = "cyc",
               preg     = "preg",
               sex      = "sex",
               all_base = c("id", "age", "bmi"),
               all_cyc  = c("id", "cyc", "preg", "opk"))




# target output ----------------------------------------------------------------

block_nm = c("beg_day_idx", "n_days", "beg_w_idx", "beg_sex_idx", "n_sex_days", "u_col", "subj_idx")

miss_age <- list(setNames(c( 0L,  5L,  0L, -1L,  0L,  1L,  0L), block_nm),
                 setNames(c( 5L,  3L,  5L,  0L,  1L,  1L,  1L), block_nm),
                 setNames(c( 8L,  7L,  8L,  1L,  2L,  1L,  2L), block_nm))

miss_bmi <- list(setNames(c(17L,  4L, 15L,  3L,  2L,  2L,  4L), block_nm),
                 setNames(c(21L,  8L, 19L,  5L,  1L,  2L,  5L), block_nm))

miss_opk <- list(setNames(c( 0L,  2L,  0L, -1L,  0L,  7L,  0L), block_nm),
                 setNames(c( 2L,  3L,  2L, -1L,  0L,  7L,  0L), block_nm),
                 setNames(c( 5L,  3L,  5L,  0L,  1L,  7L,  1L), block_nm),
                 setNames(c(13L,  2L, 13L,  2L,  1L,  7L,  2L), block_nm),
                 setNames(c(19L,  2L, 17L,  3L,  2L,  7L,  4L), block_nm))

miss_drnk <- list(setNames(c( 8L,  1L,  8L, -1L,  0L,  8L,  2L), block_nm),
                  setNames(c(11L,  1L, 11L, -1L,  0L,  8L,  2L), block_nm),
                  setNames(c(13L,  1L, 13L,  2L,  1L,  8L,  2L), block_nm),
                  setNames(c(17L,  1L, 15L, -1L,  0L,  8L,  4L), block_nm),
                  setNames(c(19L,  1L, 17L,  3L,  1L,  8L,  4L), block_nm),
                  setNames(c(25L,  1L, 23L, -1L,  0L,  8L,  5L), block_nm),
                  setNames(c(26L,  1L, 24L, -1L,  0L,  8L,  5L), block_nm),
                  setNames(c(28L,  1L, 26L, -1L,  0L,  8L,  5L), block_nm))

preg_idx <- c(rep(-1, 2), 0, 1, 2, rep(-1, 8), 3, 4, rep(-1, 10), 7, 8)
sex_miss_idx <- c(0, 1, 2, 4, 5, 6)

target_with_intercept <- list(cov_row_miss_list = list(miss_age, miss_bmi, miss_opk, miss_drnk),
                              cov_miss_w_idx    = preg_idx,
                              cov_miss_x_idx    = sex_miss_idx)

target_no_intercept <- target_with_intercept
target_no_intercept$cov_row_miss_list[[1]] <- lapply(miss_age, function(x) replace(x, 6, 0L))
target_no_intercept$cov_row_miss_list[[2]] <- lapply(miss_bmi, function(x) replace(x, 6, 1L))




# begin testing ----------------------------------------------------------------

test_get_cov_row_miss_info <- function() {

    out_with_intercept <- get_cov_row_miss_info(comb_dat,
                                                var_nm,
                                                expan_with_intercept,
                                                cov_with_intercept_info)
    checkIdentical(target_with_intercept, out_with_intercept)

    out_no_intercept <- get_cov_row_miss_info(comb_dat,
                                              var_nm,
                                              expan_no_intercept,
                                              cov_no_intercept_info)
    checkIdentical(target_no_intercept, out_no_intercept)
}
