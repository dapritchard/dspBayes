# construct data ---------------------------------------------------------------

fw_incl <- 5:8

daily <- data.frame(id  = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4),
                    cyc = c(3, 3, 4, 4, 4, 4, 4, 1, 1, 2, 2, 4, 4, 4),
                    day = c(7, 5, 1, 4, 9, 6, 8, 4, 6, 5, 8, 1, 5, 8),
                    int = rep(c("no", "yes"), 7))

cycle <- data.frame(id   = c(1, 1, 2, 2, 3, 4),
                    cyc  = c(3, 4, 1, 2, 1, 4),
                    preg = rep(c("no", "yes"), 3))

baseline <- data.frame(id  = c( 1,  3,  4,  6,  9),
                       age = c(32, 28, 26, 39, 38))

vars_nm <- list(all  = c("id", "cyc", "day", "int", "preg", "age"),
                id   = "id",
                cyc  = "cyc",
                fw   = "day",
                preg = "preg")



merge_dsp_data(daily, cycle, baseline, vars_nm, 0)
