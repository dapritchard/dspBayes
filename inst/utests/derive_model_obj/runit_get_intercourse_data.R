# construct data ---------------------------------------------------------------

# idx:     0    1    2    3    4    5    6    7    8    9   10   11   12
# -----------------------------------------------------------------------
id   <- c("a", "a", "a", "b", "b", "b", "c", "d", "d", "d", "d", "d", "e")
cyc  <- c(  1,   2,   2,   5,   5,   5,   3,   1,   1,   1,   2,   2,   7)
day  <- c(  0,   1,   2,   0,   1,   2,   1,   0,   1,   2,   0,   2,   1)
sex  <- c("y",  NA,  NA, "y",  NA, "n", "y", "n", "y", "y",  NA, "n",  NA)
yest <- c( NA, "n",  NA,  NA, "y",  NA, "n",  NA, "n", "y",  NA, "n", "n")
preg <- c("n", "y", "y", "n", "n", "n", "y", "y", "n", "n", "y", "y", "n")

# sex no missing versions.  Note that even if there are missing on the day
# before that we shouldn't consider there to be any missing in the output
sex2 <- c("y", "y", "n", "y", "n", "n", "y", "n", "y", "y", "y", "n", "y")
yes2 <- c( NA, "n", "y",  NA, "y", "n", "n",  NA, "n", "y",  NA, "n", "n")

daily <- data.frame(id, cyc, day, sex, sex_yester = yest, preg,
                    stringsAsFactors = FALSE)

daily_nomiss <- data.frame(id, cyc, day, sex = sex2, sex_yester = yes2, preg = preg,
                           stringsAsFactors = FALSE)

var_nm <- list(id = "id", cyc = "cyc", fw = "day", sex = "sex", preg = "preg")
fw_incl <- c(0, 1, 2)




# target output ----------------------------------------------------------------

# first missing in a cycle is set to 1, the remaining missing to 0
X   <- c(1L,  1L,  0L,  1L,  1L,  0L,  1L,  0L,  1L,  1L,  1L,  0L,  1L)

# note that we're using 0-based indexing here for everything
miss_cyc_1 <- c(beg_idx = 0L, n_days = 2L, subj_idx = 0L, preg_idx =  0L)
miss_cyc_2 <- c(beg_idx = 2L, n_days = 1L, subj_idx = 1L, preg_idx = -1L)
miss_cyc_5 <- c(beg_idx = 3L, n_days = 1L, subj_idx = 3L, preg_idx =  4L)
miss_cyc_6 <- c(beg_idx = 4L, n_days = 1L, subj_idx = 4L, preg_idx = -1L)
target_miss_cyc <- list(miss_cyc_1, miss_cyc_2, miss_cyc_5, miss_cyc_6)

# day_1  <- c(idx = 1L,  prev =  0L)
# day_2  <- c(idx = 2L,  prev = -1L)
# day_4  <- c(idx = 4L,  prev =  1L)
# day_10 <- c(idx = 10L, prev = -2L)
# day_12 <- c(idx = 12L, prev =  0L)
# target_miss_day <- list(day_1, day_2, day_4, day_10, day_12)
day_1  <- c(idx = 1L,  prev = 0L)
day_2  <- c(idx = 2L,  prev = 2L)
day_4  <- c(idx = 4L,  prev = 1L)
day_10 <- c(idx = 10L, prev = 2L)
day_12 <- c(idx = 12L, prev = 0L)
target_miss_day <- list(day_1, day_2, day_4, day_10, day_12)

target_sex2 <- ifelse(sex2 == "y", 1L, 0L)




# begin testing ----------------------------------------------------------------

test_no_missing <- function() {
    out <- get_intercourse_data(daily_nomiss, var_nm, fw_incl)
    checkIdentical(target_sex2, out$X)
    checkIdentical(vector("list", 0L), out$miss_cyc)
    checkIdentical(vector("list", 0L), out$miss_day)
}


test_with_missing <- function() {
    out <- get_intercourse_data(daily, var_nm, fw_incl)
    checkIdentical(X, out$X)
    checkIdentical(target_miss_cyc, out$miss_cyc)
    checkIdentical(target_miss_day, out$miss_day)
}
