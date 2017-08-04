# construct data ---------------------------------------------------------------

# create an integer, floating, character, and factor vector with which to build data
# frames
ints <- 1:6
floats <- seq(1.5, 6.5, 1)
chars <- c("a", "b", "a", "c", "a", "b")
facts <- c("one", "two", "three", "two", "three", "two") %>% factor
resp <- seq(-1, 1, length.out = 6)

orig_df <- data.frame(ints   = ints,
                      floats = floats,
                      chars  = chars,
                      facts  = facts,
                      resp   = resp,
                      stringsAsFactors = FALSE)

base_model <- formula(resp ~ ints + floats + chars + facts)
noint_model <- formula(resp ~ 0 + ints + floats + chars + facts)

expanded_base <- model.matrix(base_model, orig_df)
expanded_noint <- model.matrix(noint_model, orig_df)

get_missing_var_info(expanded_base, base_model)
get_missing_var_info(expanded_noint, noint_model)





orig_df_w_miss <- orig_df
orig_df_w_miss$ints[c(2L, 4L)] <- NA
orig_df_w_miss$chars[5L] <- NA
orig_df_w_miss$facts[c(2L, 6L)] <- NA

# option 1
current.na.action <- options("na.action")
options(na.action='na.pass')
expanded_noint_w_miss <- model.matrix(noint_model, orig_df_w_miss)
options('na.action' = current.na.action$na.action)

# option 2
model.matrix(noint_model, model.frame(noint_model, orig_df_w_miss, na.action = na.pass))

get_missing_var_info(orig_df_w_miss, noint_model)
