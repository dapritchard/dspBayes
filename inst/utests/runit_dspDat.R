load("data/RawData.RData")

id_name <- "subjId"
cyc_name <- "cycle"
sex_name <- "intercourse"
fw_name <- "reverseDay"
fw_incl <- -17:-13
model_formula <- formula(pregInd ~ 0 + factor(reverseDay) + age + bmi + gravid + cycleLen + lube)

dspDat_obj <- dspDat(model_formula,
                     baseline,
                     cycle,
                     daily,
                     id_name,
                     cyc_name,
                     sex_name,
                     fw_name,
                     fw_incl,
                     use_na       = "all",
                     req_min_days = 0L,
                     return_comb  = FALSE)
