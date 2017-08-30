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

# TODO: remove data duplication above




# target output ----------------------------------------------------------------

u_with_intercept <- expan_with_intercept
u_with_intercept[is.na(u_with_intercept[, "age"]), "age"] <- mean(u_with_intercept[, "age"], na.rm = TRUE)
u_with_intercept[is.na(u_with_intercept[, "bmi1"]), "bmi1"] <- 1
u_with_intercept[is.na(u_with_intercept[, "bmi2"]), "bmi2"] <- 0
u_with_intercept[is.na(u_with_intercept[, "bmi3"]), "bmi3"] <- 0
u_with_intercept[is.na(u_with_intercept[, "opk1"]), "opk1"] <- 1
u_with_intercept[is.na(u_with_intercept[, "drnk1"]), "drnk1"] <- 1

u_no_intercept <- expan_no_intercept
u_no_intercept[is.na(u_no_intercept[, "age"]), "age"] <- mean(u_no_intercept[, "age"], na.rm = TRUE)
u_no_intercept[is.na(u_no_intercept[, "bmi0"]), "bmi0"] <- 1
u_no_intercept[is.na(u_no_intercept[, "bmi1"]), "bmi1"] <- 0
u_no_intercept[is.na(u_no_intercept[, "bmi2"]), "bmi2"] <- 0
u_no_intercept[is.na(u_no_intercept[, "bmi3"]), "bmi3"] <- 0
u_no_intercept[is.na(u_no_intercept[, "opk1"]), "opk1"] <- 1
u_no_intercept[is.na(u_no_intercept[, "drnk1"]), "drnk1"] <- 1




# begin testing ----------------------------------------------------------------

test_get_u_miss_filled_in <- function() {

    out_with_intercept <- get_u_miss_filled_in(expan_with_intercept, cov_with_intercept_info)
    identical(u_with_intercept, out_with_intercept)

    out_no_intercept <- get_u_miss_filled_in(expan_no_intercept, cov_no_intercept_info)
    identical(u_no_intercept, out_no_intercept)
}
