# construct data ---------------------------------------------------------------

id <-  c("a", "a", "a", "b", "b", "b", "c", "d", "d", "d", "d", "d", "e")
cyc <- c(  1,   2,   2,   5,   5,   5,   3,   1,   1,   1,   2,   2,   7)

daily <- data.frame(id, cyc, stringsAsFactors = FALSE)
var_nm <- list(id = "id", cyc = "cyc")

target <- list(1L,
               2:3,
               4:6,
               7L,
               8:10,
               11:12,
               13L)


# begin testing ----------------------------------------------------------------

test_get_cyc_idx <- function() {
    checkIdentical(target, get_cyc_idx(daily, var_nm))
}
