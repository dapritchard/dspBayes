# construct data ---------------------------------------------------------------

# create an integer, character, and factor vector with which to build data
# frames
ints <- 1:6
chars <- c("a", "b", "c", "d", "e", "f")
facts <- c("one", "two", "three", "four", "five", "six") %>% factor

# create data frames with 6 rows and either 1, 2, or 3 columns
one_col_df <- data.frame(ints)
two_cols_df <- data.frame(ints, chars, stringsAsFactors = FALSE)
three_cols_df <- data.frame(ints, chars, facts, stringsAsFactors = FALSE)

# partition data frames into three two-row data frames
one_A <- one_col_df[1:2, , drop = FALSE]
one_B <- one_col_df[3:4, , drop = FALSE]
one_C <- one_col_df[5:6, , drop = FALSE]
two_A <- two_cols_df[1:2, ]
two_B <- two_cols_df[3:4, ]
two_C <- two_cols_df[5:6, ]
three_A <- three_cols_df[1:2, ]
three_B <- three_cols_df[3:4, ]
three_C <- three_cols_df[5:6, ]

# store the partitions in a list
list_of_dfs_one <- list(A=one_A, B=one_B, C=one_C)
list_of_dfs_two <- list(A=two_A, B=two_B, C=two_C)
list_of_dfs_three <- list(A=three_A, B=three_B, C=three_C)




# begin testing ----------------------------------------------------------------

# test functioning construction of data frame with one, two, or three columns
test_whole_df <- function() {
    checkEquals(rbind_similar_dfs(list_of_dfs_one), one_col_df)
    checkEquals(rbind_similar_dfs(list_of_dfs_two), two_cols_df)
    checkEquals(rbind_similar_dfs(list_of_dfs_three), three_cols_df)
}


# test length-1 input
test_length_one <- function() {
    list_of_one <- list_of_dfs_three[1L]
    checkEquals(rbind_similar_dfs(list_of_one), three_A)
}


# test length-0 input
test_no_cycles <- function() {
    # require more nonmissing days than are in the fertile window
    checkException(rbind_similar_dfs(list()), "no cycles passed the criteria")
}
