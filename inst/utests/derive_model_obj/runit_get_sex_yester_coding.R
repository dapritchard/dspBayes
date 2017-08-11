# construct data ---------------------------------------------------------------

day        <- c(  0,   1,   2,   0,   1,   2,   1,   0,   1,   2,   0,   2,   1)
sex_yester <- c( NA, "n",  NA,  NA, "y",  NA,  NA,  NA, "n", "y",  NA,  NA,  NA)
target     <- c(-2L,  0L, -1L, -2L,  1L, -1L, -1L, -2L,  0L,  1L, -2L, -1L, -1L)

daily <- data.frame(day, sex_yester)

var_nm <- list(fw = "day")
fw_incl <- c(0, 1, 2)




# begin testing ----------------------------------------------------------------

test_get_sex_yester_coding <- function() {
    checkIdentical(target, get_sex_yester_coding(daily, var_nm, fw_incl))
}
