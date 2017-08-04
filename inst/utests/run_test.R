
suite_rbind_similar_dfs <- defineTestSuite(
    name = "rbind_similar_dfs",
    dirs = file.path(getwd(), "utest"),
    testFileRegexp = "^runit_rbind_similar_dfs.R$",
    testFuncRegexp = "^test_.+")

suite_construct_template_df <- defineTestSuite(
    name = "construct_template_df",
    dirs = file.path(getwd(), "utest"),
    testFileRegexp = "^runit_construct_template_df.R$",
    testFuncRegexp = "^test_.+")

test_result <- runTestSuite(list(suite_rbind_similar_dfs,
                                 suite_construct_template_df))


printTextProtocol(test_result)
