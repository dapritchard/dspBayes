
# Partition the model variables by dataset -------------------------------------

getVarNames <- function(formula, baseline, cycle, daily, id_name, cyc_name, sex_name, fw_name) {
  explanNames <- all.vars(formula)
  pregName <- explanNames[1]
  varInSetBool <- function(set) explanNames %in% names(set)

  varNames <- list( id = id_name,
                    cyc = cyc_name,
                    preg = pregName,
                    sex = sex_name,
                    fw = fw_name,
                    basIncl = c(id_name, explanNames[varInSetBool(baseline)]),
                    cycIncl = c(id_name, cyc_name, explanNames[varInSetBool(cycle)]),
                    dayIncl = c(id_name, cyc_name, sex_name, fw_name, explanNames[varInSetBool(daily)]),
                    modelVars = list( bas = explanNames[varInSetBool(baseline)],
                                      cyc = setdiff(explanNames[varInSetBool(cycle)], pregName),
                                      day = setdiff(explanNames[varInSetBool(daily)], pregName) ) )
  return (varNames)
}




# Obtain vector of common id's that are needed ---------------------------------
#
# NULL[[id_name]] == NULL    <-- behavior of NULL object for when bas / cyc are null
# unique(NULL) == NULL
# length(NULL) == 0

getCommonId <- function(cleanDat, id_name) {
  basId <- cleanDat$bas[[id_name]]
  cycId <- unique( cleanDat$cyc[[id_name]] )
  dayId <- unique( cleanDat$day[[id_name]] )

  intrAllowNull <- function(a, b) if (!is.null(b)) intersect(a, b) else a
  idVec <- Reduce(f=intrAllowNull, x=list(dayId, cycId, basId))

  return (idVec)
}




# Obtain list of common cycles that are needed ---------------------------------

getCommonCyc <- function(cleanDat, varNames, idVec) {
  cycId <- cleanDat$cyc[[varNames$id]]
  dayId <- cleanDat$day[[varNames$id]]
  cyc_name <- varNames$cyc
  n <- length(idVec)

  if (is.null(cleanDat$cyc))
    cycList <- lapply(idVec, function(x) unique(cleanDat$day[x == dayId, cyc_name]))
  else {
    cycById <- list( cyc = lapply(idVec, function(x) cleanDat$cyc[x == cycId, cyc_name]),
                     day = lapply(idVec, function(x) cleanDat$day[x == dayId, cyc_name]) )
    cycList <- lapply(1:n, function(j) intersect(cycById$cyc[[j]], cycById$day[[j]]) )
  }
  return (cycList)
}




# Convert possible 0/1 or "yes"/"no" to boolean --------------------------------
#
# PRE: already checked that all values are logical or numeric or start with "y","Y","n", or "N"

convToBool <- function(x, keepNA=TRUE) {
  convChar <- function(c) {
    if (is.na(c)) NA
    else if (identical(c, "n") || identical(c, "N")) FALSE
    else TRUE
  }

  if (is.logical(x)) boolVec <- x
  else if (is.numeric(x)) boolVec <- as.logical(x)
  else boolVec <- sapply(substr(x, 1, 1), convChar)

  return ( replace(boolVec, is.na(boolVec), ifelse(keepNA, TRUE, FALSE)) )
}
