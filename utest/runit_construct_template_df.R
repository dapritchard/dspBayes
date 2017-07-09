# construct data ---------------------------------------------------------------

# create an integer, character, and factor vector with which to build data
# frames
ints <- 1:6
chars <- c("a", "b", "c", "d", "e", "f")
facts <- c("one", "two", "three", "four", "five", "six") %>% factor

# all NA versions of `ints`, `chars`, and `facts`
ints_NA <- rep(NA_integer_, 6L)
chars_NA <- rep(NA_character_, 6L)
facts_NA <- rep(NA_integer_, 6L)
attributes(facts_NA) <- attributes(facts)

# create data frames with 6 rows and either 1, 2, or 3 columns
one_col_df <- data.frame(ints)
two_cols_df <- data.frame(ints, chars, stringsAsFactors = FALSE)
three_cols_df <- data.frame(ints, chars, facts, stringsAsFactors = FALSE)

# all NA versions of `one_col_df`, `two_col_df`, `three_col_df`
one_col_NA <- data.frame(ints_NA) %>% setNames(., names(one_col_df))
names(one_col_NA) <- names(one_col_df)
two_cols_NA <- data.frame(ints_NA, chars_NA, stringsAsFactors = FALSE)
names(two_cols_NA) <- names(two_cols_df)
three_cols_NA <- data.frame(ints_NA, chars_NA, facts_NA, stringsAsFactors = FALSE)
names(three_cols_NA) <- names(three_cols_df)




# construct data ---------------------------------------------------------------

test_template <- function() {
    checkEquals(construct_template_df(one_col_df, 1:6), one_col_NA)
    checkEquals(construct_template_df(two_cols_df, 1:6), two_cols_NA)
    checkEquals(construct_template_df(three_cols_df, 1:6), three_cols_NA)
}


test_various_row_sizes <- function() {

    # only one row and one column
    checkEquals(construct_template_df(one_col_df[1L, , drop = FALSE], 1L),
                one_col_NA[1L, , drop = FALSE])

    # second argument has smaller length then number of rows in first argument
    checkEquals(construct_template_df(three_cols_df, 1:3), three_cols_NA[1:3, ])

    # second argument has larger length then number of rows in first argument
    checkEquals(construct_template_df(three_cols_df[1:3, ], 1:6), three_cols_NA)
}
