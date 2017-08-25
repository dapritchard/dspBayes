# TODO: test for tracking variables (late addition)

# construct data ---------------------------------------------------------------

# store as variables to ensure consistency between baseline, cycle, and daily
# datasets
id <- c(1, 2)
cyc <- c(1, 2, 1, 2)
preg_status <- c("y", "n", "y", "y")

# basic data variables
id_name <- "id"
cyc_name <- "cyc"
sex_name <- "sex"
fw_name <- "day"
preg_name <- "preg"
fw_days <- c(2, 3, 4)

# construct daily data.  `daily` is without a pregnancy variable, and
# `daily_preg` has one added.
daily <- data.frame(id       = rep(id, each = 10),
                    cyc      = rep(cyc, each = 5),
                    day      = rep(1:5, 4),
                    sex      = rep(c("no", "yes"), 5) %>% factor,
                    exercise = rep(c("none", "cardio", "weights"), length.out = 10))
daily_preg <- data.frame(daily, preg = rep(preg_status, each = 5))

# construct cycle data
cycle <- data.frame(id      = rep(id, each = 2),
                    cyc     = cyc,
                    preg    = preg_status,
                    opk_use = c("No", NA, "Yes", "No"))

# construct baseline data
baseline <- data.frame(id  = id,
                       age = factor(c("<35", "38+"), levels = c("<35", "35-38", "38+")),
                       bmi = c(24.4, 20.2))

# construct 3 different models of complexity and corresponding `var_nm` objects
extract_wrapper <- function(input_model) {
    extract_var_nm(input_model,
                   baseline,
                   cycle,
                   daily,
                   id_name,
                   cyc_name,
                   sex_name,
                   fw_name)
}
var_nm_minimal <- extract_wrapper(preg ~ 0 + day)
var_nm_age <- extract_wrapper(preg ~ 0 + day + age)
var_nm_all <- extract_wrapper(preg ~ 0 + day + exercise + opk_use + age + bmi)

# arbitrary permutations of the data
baseline_perm <- baseline[c(2, 1), ]
cycle_perm <- cycle[c(3, 1, 4, 2), ]
day_perm_idx <- c(14, 3, 17, 11, 6, 8, 1, 19, 5, 13, 2, 16, 7, 9, 10, 20, 15, 18, 4, 12)
daily_perm <- daily[day_perm_idx, ]

# remove some daily observations.  Note that an entire cycle is missing, which
# should be recovered if we have a corresponding observation in the cycle data.
day_nonexist_idx <- c(1:5, 11, 13, 15, 17, 19)
daily_nonexist <- daily[setdiff(day_perm_idx, day_nonexist_idx), ]




# testing objects --------------------------------------------------------------

# pick off days 2-4 out of each 5-day block
day_keep_idx <- c(2:4, 7:9, 12:14, 17:19)

# each cycle has 3 fertile window days
cyc_keep_idx <- rep(1:4, each = 3)

# each subject has 2 cycles of 3 fertile window days each
base_keep_idx <- rep(1:2, each = 6)

# target daily data, in the order in which they appear
target_id       <- daily[day_keep_idx, id_name]
target_cyc      <- daily[day_keep_idx, cyc_name]
target_day      <- daily[day_keep_idx, fw_name]
target_sex      <- daily[day_keep_idx, sex_name]
target_exercise <- daily[day_keep_idx, "exercise"]

# target cycle data, in the order in which they appear
target_preg    <- cycle[cyc_keep_idx, preg_name]
target_opk_use <- cycle[cyc_keep_idx, "opk_use"]

# target baseline data, in the order in which they appear
target_age <- baseline[base_keep_idx, "age"] %>% droplevels
target_bmi <- baseline[base_keep_idx, "bmi"]

# target merged data using all covariates.  They have to be ordered first by
# dataset (daily, cycle, baseline), and then by the order of the columns within
# the dataset.
target_all <- data.frame(id       = target_id,
                         cyc      = target_cyc,
                         day      = target_day,
                         sex      = target_sex,
                         exercise = target_exercise,
                         preg     = target_preg,
                         opk_use  = target_opk_use,
                         age      = target_age,
                         bmi      = target_bmi,
                         stringsAsFactors = FALSE)

# column indices corresponding to `var_nm_minimal` and `var_nm_age` models
extract_col_idx <- function(nm, target = target_all) {
    sapply(nm, function(x) which(x == colnames(target)))
}
min_col_idx <- extract_col_idx(c(id_name, cyc_name, fw_name, sex_name, preg_name))
age_col_idx <- extract_col_idx(c(id_name, cyc_name, fw_name, sex_name, preg_name, "age"))

# target merged data corresponding to `var_nm_minimal` and `var_nm_age` models
target_minimal <- target_all[, min_col_idx]
target_age <- target_all[, age_col_idx]

# target merged data using all covariates when days are nonexistent in the
# fertile window
target_nonexist_idx <- which(day_keep_idx %in% day_nonexist_idx)
target_nonexist <- target_all
target_nonexist[target_nonexist_idx, c(sex_name, "exercise")] <- NA

# create daily data with missing.  By setting the missing to the same days
# as we used for the nonexistent data we get the same results, since the
# function throws out those days.
daily_w_id_miss <- daily_w_cyc_miss <- daily_w_fw_miss <- daily
daily_w_id_miss[day_nonexist_idx, id_name] <- NA
daily_w_cyc_miss[day_nonexist_idx, cyc_name] <- NA
daily_w_fw_miss[day_nonexist_idx, fw_name] <- NA

# create subject data with missing.  Both of these correspond to throwing
# out observations 1-3 of the `target_all` data.
cycle_w_id_miss <- cycle_w_cyc_miss <- cycle_w_preg_miss <- cycle
cycle_w_id_miss[1L, id_name] <- NA
cycle_w_cyc_miss[1L, cyc_name] <- NA
cycle_w_preg_miss[1L, preg_name] <- NA
target_miss_cyc_1 <- target_all[-c(1:3), ] %>% structure(., row.names = 1:9)




# `merge_dsp_data` function wrappers -------------------------------------------

merge_defaults <- function(baseline     = baseline_perm,
                           cycle        = cycle_perm,
                           daily        = daily_perm,
                           var_nm       = var_nm_minimal,
                           fw_incl      = fw_days,
                           min_days_req = 0L) {
    merge_dsp_data(baseline, cycle, daily, var_nm, fw_incl, min_days_req)
}


# for testing against various choices of `var_nm`.  Uses the permuted data as
# arguments.
merge_model <- function(var_nm,
                        baseline     = baseline_perm,
                        cycle        = cycle_perm,
                        daily        = daily_perm,
                        fw_incl      = fw_days,
                        min_days_req = 0L) {
    merge_dsp_data(baseline, cycle, daily, var_nm, fw_incl, min_days_req)
}


# for testing when daily has nonexistent entries for some of the days in the
# fertile window, with various choices of `min_days_req`
merge_nonexist_days <- function(min_days_req,
                                baseline     = baseline_perm,
                                cycle        = cycle_perm,
                                daily        = daily_nonexist,
                                var_nm       = var_nm_all,
                                fw_incl      = fw_days) {
    merge_dsp_data(baseline, cycle, daily, var_nm, fw_incl, min_days_req)
}




# begin testing ----------------------------------------------------------------

# various versions of `var_nm` are determined through calls to `extract_var_nm`
# which depend on the model specification
test_model_specifications <- function() {

    out_minimal <- merge_model(var_nm_minimal)
    checkEquals(target_minimal, out_minimal)

    out_age <- merge_model(var_nm_age)
    checkEquals(target_age, out_age)

    out_all <- merge_model(var_nm_all)
    checkEquals(target_all, out_all)
}


# when days in the fertile window are nonexistent in the daily data then we want
# to fill them in and insert an NA for the corresponding data
test_nonexistent_fw_days <- function() {

    out_min_0 <- merge_nonexist_days(0)
    checkEquals(target_nonexist, out_min_0)

    out_min_1 <- merge_nonexist_days(1)
    target <- target_nonexist[-c(1:3), ] %>% structure(., row.names = 1:9)
    checkEquals(target, out_min_1)

    out_min_2 <- merge_nonexist_days(2)
    target <- target_nonexist[-c(1:3, 10:12), ] %>% structure(., row.names = 1:6)
    target$opk_use <- droplevels(target$opk_use)
    checkEquals(target, out_min_2)

    checkException(merge_nonexist_days(99))
}


# should work the same if pregnancy data is in daily rather than cycle
test_preg_in_daily <- function() {

    cyc_preg_idx <- which(colnames(cycle_perm) == preg_name)
    out <- merge_model(var_nm_all, daily = daily_preg, cycle = cycle[, -cyc_preg_idx])
    checkEquals(target_all, out)
}


# if we have missing in any of the subject id, cycle number, pregnancy, or
# fertile window value then we should throw out this data
test_missing_id_cyc_fw <- function() {

    out <- merge_model(var_nm_all, daily = daily_w_id_miss)
    checkEquals(target_nonexist, out)

    out <- merge_model(var_nm_all, daily = daily_w_cyc_miss)
    checkEquals(target_nonexist, out)

    out <- merge_model(var_nm_all, daily = daily_w_fw_miss)
    checkEquals(target_nonexist, out)

    out <- merge_model(var_nm_all, cycle = cycle_w_id_miss)
    checkEquals(target_miss_cyc_1, out)

    out <- merge_model(var_nm_all, cycle = cycle_w_cyc_miss)
    checkEquals(target_miss_cyc_1, out)

    out <- merge_model(var_nm_all, cycle = cycle_w_preg_miss)
    checkEquals(target_miss_cyc_1, out)
}


# baseline and cycle datasets are allowed to be NULL
test_null_dataset <- function() {

    out <- merge_defaults(baseline = NULL, cycle = NULL, daily = daily_preg)
    checkEquals(target_minimal, out)

    out <- merge_defaults(cycle = NULL, daily = daily_preg, var_nm = var_nm_age)
    checkEquals(target_age, out)
}
