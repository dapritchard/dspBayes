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




# target output ----------------------------------------------------------------

get_means <- function(expan_df, col_idx) {

    categ_means <- apply(expan_df[, col_idx, drop = FALSE], 2, mean, na.rm = TRUE) %>% setNames(., NULL)

    if (sum(categ_means) == 1) {
        return(c(categ_means))
    } else {
        return(c(categ_means, 1 - sum(categ_means)))
    }
}

target_with_intercept <- list(age  = list(idx             = 2L,
                                          composition_nm  = "age",
                                          categ           = FALSE,
                                          empirical_probs = numeric(0L),
                                          n_categs        = 1L),
                              bmi  = list(idx             = 3:5,
                                          composition_nm  = "bmi",
                                          categ           = TRUE,
                                          empirical_probs = get_means(expan_with_intercept, 3:5),
                                          n_categs        = 4L),
                              opk  = list(idx             = 8L,
                                          composition_nm  = "opk",
                                          categ           = TRUE,
                                          empirical_probs = get_means(expan_with_intercept, 8L),
                                          n_categs        = 2L),
                              drnk = list(idx             = 9L,
                                          composition_nm  = "drnk",
                                          categ           = TRUE,
                                          empirical_probs = get_means(expan_with_intercept, 9L),
                                          n_categs        = 2L))

target_no_intercept <- target_with_intercept
target_no_intercept$age$idx <- 1L
target_no_intercept$bmi$idx <- 2:5
target_no_intercept$bmi$empirical_probs <- get_means(expan_no_intercept, 2:5)




# begin testing ----------------------------------------------------------------

test_get_cov_col_miss_info <- function() {

    out_with_intercept <- get_cov_col_miss_info(expan_with_intercept , model_with_intercept, "all")
    checkIdentical(target_with_intercept, out_with_intercept)

    out_no_intercept <- get_cov_col_miss_info(expan_no_intercept , model_no_intercept, "all")
    checkIdentical(target_no_intercept, out_no_intercept)
}
