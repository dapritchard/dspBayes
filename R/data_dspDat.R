#' Specify model variables for the day-specific probabilities MCMC sampler
#'
#' \code{dspDat} is used to create an object of \code{class} \code{"dspDat"};
#' the resultant object may then be used as input to the \code{dsp} function to
#' sample an MCMC chain for the methodology proposed by Dunson and Stanford in
#' \emph{Bayesian Inferences on Predictors of Conception Probabilities} (2005).
#' The \code{dspDat} function is essentially a convenience function provided to
#' (if necessary) merge multiple datasets of varying time-specificities, as is
#' common for the type of fertility data for which the aformentioned methodology
#' is designed.
#'
#' The \code{class} \code{"dspDat"} is equipped with a \code{summary} function.
#'
#' @param dsp_model An object of class \code{\link[stats]{formula}} (or one that
#'     can be coerced to that class).  The term on the left-hand side of the
#'     formula must be the name of a column in either \code{cycle} or
#'     \code{daily}, with the observations in this column indicating whether the
#'     cycle (or cycle in which the day was a part of) resulted in a successful
#'     pregnancy.  The terms on the right-hand side of the formula must each be
#'     a name of a column in any of one of \code{baseline}, \code{cycle}, or
#'     \code{daily}.
#'
#' @param baseline Either \code{NULL} or a \code{data.frame} (or an object that
#'     can be coerced to \code{data.frame}).  If non-\code{NULL}, then contains
#'     data in the form of one observation (row) per study subject, and may
#'     include covariates (columns) such as e.g. age at time entering study, BMI
#'     at time entering study, gravidity status, etc. A value of \code{NULL} may
#'     mean that either no such variables are to be included in the model, or
#'     that such data has already been expanded and is included in the
#'     \code{cycle} or \code{daily} data.  If a non-\code{NULL} object is
#'     supplied, then \code{baseline} must include a column with name as
#'     specified by the \code{id_name} parameter which provides a study ID for
#'     each observation.
#'
#' @param cycle Either \code{NULL} or a \code{data.frame} (or an object that can
#'     be coerced to \code{data.frame}).  If non-\code{NULL}, then contains data
#'     in the form of one observation (row) per study cycle, and may include
#'     columns such as e.g. cycle pregnancy indicator, attempt cycle number, or
#'     cycle length, etc. A value of \code{NULL} indicates that such data has
#'     already been expanded and is included in the \code{daily} data.  If a
#'     non-\code{NULL} object is supplied, then \code{cycle} must include a
#'     column with name as specified by the \code{id_name} parameter which
#'     provides a study id for each observation, and must include a column with
#'     name as specified by \code{cyc_name} which provides a cycle number for
#'     each observation.  If the pregnancy outcome variable is included in this
#'     data, then the column must have the name as specified by the left-hand
#'     term in the \code{dsp_model} parameter.
#'
#' @param daily A \code{data.frame} (or an object that can be coerced to
#'     \code{data.frame}).  Contains data in the form of one observation (row)
#'     per study day.  May include data such as e.g. day intermenstual bleeding
#'     indicator, day cervical mucus type etc.  Must include a column with name
#'     as specified by the \code{id_name} parameter which provides a study id
#'     for each observation, a column with name as specified by \code{cyc_name}
#'     which provides a cycle number for each observation, and a column with
#'     name as specified by \code{sex_name} which provides an indicator of
#'     whether the observation (day) is within the fertile window.  If a cycle
#'     pregnancy outcome column was not provided in the \code{cycle} data, then
#'     one must be provided in the \code{daily} data, and must have name as
#'     specified by the left-hand term in the \code{dsp_model} parameter.
#'
#' @param id_name A string specifying the name of the column in each of the
#'     non-\code{NULL} \code{baseline}, \code{cycle}, or \code{daily} objects
#'     such that the column observations provide the study id for the subject to
#'     which each observation belongs to.  The name of the column must be the
#'     same for each of the non-\code{NULL} datasets.
#'
#' @param cyc_name A string specifying the name of the column in each of the
#'     non-\code{NULL} \code{cycle} or \code{daily} objects such that the column
#'     observations provide the cycle number to which each observation belongs
#'     to.  If \code{cycle} is non-\code{NULL}, then the name of the column must
#'     be the same for both the \code{cycle} and \code{daily} data.
#'
#' @param sex_name A string specifying the name of the column in the
#'     \code{daily} data such that the column observations provide an indicator
#'     of whether unprotected vaginal intercourse occurred during that day.
#'
#' @param fw_name If non-\code{NULL}, then a string specifying the name of the
#'     column in the \code{daily} data such that the column observations provide
#'     an indicator of whether the observation is part of a cycle's fertile
#'     window. If \code{NULL}, then it is assumed that the \code{daily} data has
#'     already been restricted to only observations that occurred during each
#'     cycle's fertile window.
#'
#'     As a convenience, if the name specified by \code{fw_name} is included in
#'     the \code{dsp_model} parameter, then it is interpreted to mean that a
#'     factor variable is to be included in the model corresponding to the day
#'     number in the fertile window, i.e. fertile window day 1, fertile window
#'     day 2, etc.  Warning: this assumes that the observations within a cycle
#'     are in chronological order.  TODO: needs updated
#'
#' @param fw_incl An atomic vector specifying the days in the cycle that are in
#'     the fertile window.  In more detail, the elements of \code{fw_incl}
#'     should specify which elements of the cycle day variable (as specified by
#'     \code{fw_name}) are a part of the fertile window.  So for example, if the
#'     cycle day data is represented by numeric values \eqn{ 1, 2, \dots }, and
#'     we want to specify days 13 through 17 as belonging to the fertile window,
#'     then this can be done using e.g. the argument \code{13:17}.
#'
#' @param use_na One of either \code{"none"} or \code{"sex"}.  If \code{"none"}
#'     then observations with missing data are removed from the model.  If
#'     \code{"sex"} then observations with missing intercourse data are included
#'     in the model conditional on no other data missing in the observation.
#'     See \emph{Data Processing Steps} for more details.  TODO: update
#'
#' @param req_min_days A single nonnegative numeric value specifying a minimum
#'     number of day-specific observations in the data for a given cycle which
#'     we require to include a cycle in the model.  So for example, specifying
#'     the value 0 means that we do not require any day-specific data in order
#'     to include a cycle in the model as long as we have the pregnancy outcome
#'     status for the cycle.  Specifying a value of 1 would require at least 1
#'     day-specific observation in the data for a given cycle for the cycle to
#'     be included, and so on.
#'
#'     Note that this parameter only has an effect when \code{use_na} has a
#'     value of either \code{"intercourse"} or \code{"all"}.  For if not, then
#'     cycles without a day-specific observation for every day in the fertile
#'     window will have missing values for the intercourse data, which will
#'     cause the cycle to be removed.
#'
#' @param keep_data Either \code{TRUE} or \code{FALSE}, specifying whether the
#'     merged data should be included in the return object.  In broad terms,
#'     there are three main data processing steps performed by the \code{dspDat}
#'     function: merging the datasets and remove unneeded variables,
#'     conditionally removing cycles with missing data, and removing days from
#'     the data in which intercourse did not occur.  If \code{TRUE} is
#'     specified, then the return object will include elements named
#'     \code{comb_dat}, \code{clean_dat}, and \code{sex_days_dat} which are the
#'     datasets for the various steps in the process, provided as
#'     \code{data.frame}s.
#'
#'     It is strongly recommended that the first time a model is specified using
#'     \code{dspDat}, that the user visually inspect these datasets to ensure
#'     that the data processing steps are being performed as expected.
#'
#' @details It is natural to record fertility study data in up to three datasets
#'     of varying time-specificities.  First, a dataset of variables that do not
#'     change throughout the study which we denote as the \code{baseline} data,
#'     second a dataset of cycle-specific variables which we denote as the
#'     \code{cycle} data, and third a dataset of day-specific variables which we
#'     denote as the \code{daily} data.  \code{dspDat} is provided as a
#'     convenience function which merges all of the provided datasets into one
#'     day-specific dataset and creates some internal objects for use by the
#'     MCMC sampler function \code{dsp}.
#'
#'     At a minimum the \code{daily} data must be provided so that daily
#'     intercourse data is available. \code{baseline} and \code{cycle} data are
#'     optional, so long as pregnancy information is included in one of either
#'     the \code{cycle} data or \code{daily} data.  For example, if the data was
#'     collected only in a daily format or has already been combined, then only
#'     a day-specific dataset would need to be passed to \code{dspDat}.
#'
#'     The usual \code{\link[stats]{model.matrix}} is used to construct the
#'     design matrix for the specified model, so any of the usual
#'     \code{\link[stats]{formula}} commands are available. In particular, a
#'     formula has an implied intercept term which may not be desireable for
#'     these types of models.  To remove this use either e.g. \code{y ~ x - 1}
#'     or \code{y ~ 0 + x}.
#'
#' @return \code{dspDat} returns an object of \code{\link[base]{class}}
#'     \code{"dspDat"}. An object of class \code{"dspDat"} is a list containing
#'     the following components:
#'
#'     \describe{
#'
#'         \item{TODO}{this section needs reworked}
#'
#'         \item{\code{cleanDat}}{A list containing objects \code{bas},
#'         \code{cyc}, and \code{day}, which are the datasets after removing
#'         missing and reducing the \code{daily} data to fertile window days as
#'         described in \emph{Data Processing Steps}. If \code{NULL} was
#'         supplied for \code{baseline} or \code{cycle}, then the value of
#'         \code{bas} or \code{cyc} is also \code{NULL}. }
#'
#'         \item{\code{redDat}}{A list containing objects \code{bas},
#'         \code{cyc}, and \code{day}, which are the datasets after reducing the
#'         cleaned data to the set of IDs and cycles that are common to every
#'         non-\code{NULL} dataset.  If \code{NULL} was supplied for
#'         \code{baseline} or \code{cycle}, then the value of \code{bas} or
#'         \code{cyc} is also \code{NULL}. }
#'
#'         \item{\code{combDat}}{ ******* }
#'
#'         \item{\code{modelObj}}{A list containing objects \code{Y}, \code{X},
#'         \code{U}, and \code{id}. \code{Y}, \code{X}, \code{U} are as in the
#'         Dunson and Stanford paper, and \code{id} is a vector of subject IDs
#'         such that each observation specifies the subject ID for the
#'         corresponding observation. }
#'
#'         \item{\code{samplerObj}}{A list containing objects for use by the
#'         \code{dsp} function when executing the MCMC algorithm}
#'
#'         \item{\code{datInfo}}{A list containing objects for use by the
#'         \code{\link{summary}} function}
#'
#'     }
#'
#' @section Data Processing Steps:
#'
#'     \describe{
#'
#'         \item{Cleaning data}{If either a \code{baseline} or \code{cycle}
#'         dataset is provided, then all observations that contain missing data
#'         among the model variables are removed.  All non-fertile window days
#'         are removed from the \code{daily} dataset, and any cycles that either
#'         contain missing in the fertile window or have too many or too few
#'         fertile window days are also removed.}
#'
#'         \item{Reducing data}{Each non-\code{NULL} dataset is reduced to the
#'         set of IDs and cycles that are common to every non-\code{NULL}
#'         dataset.}
#'
#'     }
#'
#' @author David Pritchard
#'
#' @references Dunson, David B., and Joseph
#'     B. Stanford. "Bayesian inferences on predictors of conception probabilities."
#'     \emph{Biometrics} 61.1 (2005): 126-133.
#'
#' @export




# Mung data into format for use by mcmc sampler ================================

dspDat <- function(dsp_model,
                   baseline      = NULL,
                   cycle         = NULL,
                   daily,
                   id_name,
                   cyc_name,
                   sex_name,
                   fw_name,
                   fw_incl,
                   fw_day_before = NULL,
                   use_na        = "none",
                   req_min_days  = 0L,
                   keep_data     = TRUE) {

    # TODO: check valid input:

    var_nm <- extract_var_nm(dsp_model,
                             baseline,
                             cycle,
                             daily,
                             id_name,
                             cyc_name,
                             sex_name,
                             fw_name)

    # merge the data provided by `baseline`, `cycle`, and `daily` into a data
    # frame
    comb_dat <- merge_dsp_data(baseline, cycle, daily, var_nm, fw_incl, fw_day_before, req_min_days)

    # conditionally remove any cycles from `comb_dat` that have either missing data
    # in the intercourse variable, other covariates, or both, depending on the value
    # of `use_na`
    clean_dat <- remove_cycs_with_miss(comb_dat, var_nm, fw_incl, use_na)

    # fit a logistic regression model with intercourse as the outcome variable
    # and sex yesterday and the model covariates as the predictor variables.
    # Note that we have to do this before performing `remove_days_no_sex` or
    # else we lose the "no intercourse" observations from the model
    tau_fit <- get_tau_fit(comb_dat, var_nm, dsp_model, use_na)

    # removes all observations from `comb_dat` for which intercourse did not occur.
    # Observations with a missing value for intercourse remain in the data.
    sex_only_dat <- remove_days_no_sex(comb_dat, var_nm)

    #### TODO check if data is collinear or constant within outcome ####

    dsp_data <- derive_model_obj(sex_only_dat, var_nm, fw_incl, dsp_model, use_na, tau_fit)

    # TODO: Stats related to munging process for use by summary fcn
    # datInfo <- getDatInfo(dsp_model, baseline, cycle, daily, cleanDat,
    #                       redDat, modelObj, varNames, fw_len, idVec, cycList)

    # conditionally add the data munging steps to the return object
    if (keep_data) {
        dsp_data$comb_dat <- comb_dat
        dsp_data$clean_dat <- clean_dat
        dsp_data$sex_only_dat <- sex_only_dat
    }

    structure(dsp_data, class="dspDat")
}




# Describes the data munging process ===========================================

summary.dspDat <- function(dspDat) {
    datInfo <- dspDat$datInfo
    hline <- paste0(rep("-", 60), collapse="")

    # internal functions to assist printing --------------------------------------

    printStats <- function(title, dat) {
        cat("\n", hline, "\n", title, ":\n\n",
            "    baseline data:  ", dat$bas$sub, " subjects\n",
            "    cycle data:     ", dat$cyc$sub, " subjects with ", dat$cyc$cyc, " cycles\n",
            "    daily data:     ", dat$day$sub, " subjects with ", dat$day$cyc, " cycles and ",
            dat$day$day, " days\n", sep="")
    }

    printVars <- function(charVec, returnWidth=40) {
        if (is.null(charVec))
            return (NULL)

        # else
        currLen <- 0
        outVec <- NULL

        for (i in 1:length(charVec)) {

            if (currLen >= returnWidth) {
                outVec <- paste0(outVec, "\n", paste(rep("", 22), collapse=" "), charVec[i], ", ")
                currLen <- 0
            }
            else {
                outVec <- paste0(outVec, charVec[i], ", ")
                currLen <- currLen + nchar(charVec[i])
            }
        }

        return ( substr(outVec, start=1, stop=(nchar(outVec) - 2)) )
    }

    # ----------------------------------------------------------------------------

    printStats("Raw data", datInfo$numRaw)
    printStats("Clean data", datInfo$numClean)

    numRed <- datInfo$numRed
    cat("\n", hline, "\nCombining the data:\n\n",
        "    common observations:  ", numRed$sub,
        " subjects with ", numRed$cyc,
        " cycles and ", numRed$day, " days\n", sep="")

    cat("\n", hline, "\nThe model variables:\n\n",
        "    variable names:  ", printVars(datInfo$modelVars), "\n",
        "    design matrix:   ", printVars(datInfo$designMatVars), "\n\n", sep="")

    # TODO: ave number of cycles in study (tot, preg, not preg), num pregnant, num sex

    cat(hline, "\n", sep="")
}




# Check if class of object is "dspDat" -----------------------------------------

is.dspDat <- function(dat) {
    identical(attributes(dspDat)$class, "dspDat")
}
