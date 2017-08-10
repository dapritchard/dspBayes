# construct data ---------------------------------------------------------------

# idx:    0    1    2    3    4    5    6    7    8    9   10   11   12
# -----------------------------------------------------------------------
id  <- c("a", "a", "a", "b", "b", "b", "c", "d", "d", "d", "d", "d", "e")
cyc <- c( 1,   2,   2,   5,   5,   5,   3,   1,   1,   1,   2,   2,   7)
day <- c(-1,   0,   1,  -1,   0,   1,   0,  -1,   0,   1,  -1,   1,   0)
sex <- c( 1,  NA,  NA,   1,  NA,   0,   1,   0,   1,   1,  NA,   0,  NA)
yst <- c(2L,  0L,  2L,  2L,  1L,  2L,  2L,  2L,  0L,  1L,  2L,  2L,  2L)

sex_nomiss <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0)

daily <- data.frame(id, cyc, day, sex, sex_yester = yst, stringsAsFactors = FALSE)
daily_nomiss <- data.frame(id, cyc, day, sex_nomiss, sex_yester = yst, stringsAsFactors = FALSE)
colnames(daily_nomiss) <- colnames(daily)

var_nm <- list(id = "id", cyc = "cyc", fw = "day", sex = "sex")




# target output ----------------------------------------------------------------

# first missing in a cycle is set to 1, the remaining missing to 0
X   <- c(1L,  1L,  0L,  1L,  1L,  0L,  1L,  0L,  1L,  1L,  1L,  0L,  1L)

# note that we're using 0-based indexing here for everything
miss_cyc_1 <- list(beg_idx =  1L, n_days = 2L, subj_idx = 0L, n_miss = 2L)
miss_cyc_2 <- list(beg_idx =  3L, n_days = 3L, subj_idx = 1L, n_miss = 1L)
miss_cyc_5 <- list(beg_idx = 10L, n_days = 2L, subj_idx = 3L, n_miss = 1L)
miss_cyc_6 <- list(beg_idx = 12L, n_days = 1L, subj_idx = 4L, n_miss = 1L)
target_miss_cyc <- list(miss_cyc_1, miss_cyc_2, miss_cyc_5, miss_cyc_6)

day_1  <- list(idx = 1L,  prev = 0L)
day_2  <- list(idx = 2L,  prev = 2L)
day_4  <- list(idx = 4L,  prev = 1L)
day_10 <- list(idx = 10L, prev = 2L)
day_12 <- list(idx = 12L, prev = 2L)
target_miss_day <- list(day_1, day_2, day_4, day_10, day_12)

n_max_miss <- 2L




# begin testing ----------------------------------------------------------------

test_no_missing <- function() {
    out <- get_intercourse_data(daily_nomiss, var_nm)
    checkIdentical(as.integer(sex_nomiss), out$X)
    checkIdentical(vector("list", 0L), out$miss_cyc)
    checkIdentical(vector("list", 0L), out$miss_day)
    checkIdentical(0L, out$n_max_miss)
}


test_with_missing <- function() {
    out <- get_intercourse_data(daily, var_nm)
    checkIdentical(X, out$X)
    checkIdentical(target_miss_cyc, out$miss_cyc)
    checkIdentical(target_miss_day, out$miss_day)
    checkIdentical(2L, out$n_max_miss)
}
