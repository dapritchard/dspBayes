# calculate the density function for the skew normal distribution using the
# definition found at https://en.wikipedia.org/wiki/Skew_normal_distribution
dskewnorm <- function(x, mu = 0, sigma = 1, shape = 0) {
    x_std <- (x - mu) / sigma
    2 / sigma * dnorm(x_std) * pnorm(shape * x_std)
}


# calculate the density function for a scaled version of the skew normal
# distribution C * f(x), such that C is chosen so that `max_x { C * f(x) } = 1`
standardized_dskewnorm <- function(x, sigma = 1, shape = 0) {
    mode_info <- find_skewnorn_mode(sigma, shape)
    # sometimes the mode can't be found
    if (is.atomic(mode_info)) {
        return(rep(0, length(x)))
    }
    normalizing_const <- mode_info$objective
    shift_const <- mode_info$maximum
    dskewnorm(x + shift_const, 0, sigma, shape) / normalizing_const
}


# find the maximum value of the skewnorm density function for the specified
# parameters
find_skewnorn_mode <- function(sigma = 1, shape = 0) {
    f <- function(x) dskewnorm(x, 0, sigma, shape)
    if (sigma <= 0) return(0)
    optimize(f, c(-10, 10), maximum = TRUE)
}


# find the parameters that minimize the squared-error loss between the observed
# values `y` and the predicted values y, where predicted y is constrained to be
# follow a standardized skew normal distribution
find_optimal_params <- function(x, y, starting = c(sigma = 1, shape = 0)) {

    # calculate the objective function for the choice of standardized skew
    # normal distribution using the parameters specified by `params`
    calc_scaled_dskew_objective <- function(params) {

        sigma = params[1L]
        shape = params[2L]

        y_hat <- standardized_dskewnorm(x, sigma = sigma, shape = shape)
        crossprod(y - y_hat)
    }

    optim(starting, calc_scaled_dskew_objective)
}


# fit the empirical curves provided in `probs_df`, for each triplet of starting
# values obtained by taking the cartesian product of all the parameter vectors
grid_search_optimal_params <- function(probs_df, grid_vals_sigma, grid_vals_shape) {

    grid_df <- expand.grid(grid_vals_sigma, grid_vals_shape)
    grid_mat <- as.matrix(grid_df)

    opt <- apply(grid_mat, 1, find_optimal_params, x = probs_df$x, y = probs_df$y)
    opt_df <- tibble::as_tibble(purrr::transpose(opt))

    tibble::tibble(
                grid_vals_sigma = grid_df[, 1L],
                grid_vals_shape = grid_df[, 2L],
                convergence     = (unlist(opt_df$convergence) == 0),
                objective       = unlist(opt_df$value),
                n_iter          = opt_df$counts,
                message         = opt_df$message,
                param_vals      = opt_df$par
            )
}




# standardize data -------------------------------------------------------------

standardize_probs <- function(df) {
    max_idx <- which.max(df$y)
    df$x <- df$x - df$x[max_idx]
    df$y <- df$y / df$y[max_idx]
    df
}


# calculate the sum of squares distances between each parameter vector in
# comparison to all of the other vectors
calc_parameter_dist <- function(params_mat) {
    sum_of_sq_dist <- function(x) {
        sum(apply(params_mat, 1L, function(y) crossprod(y - x)))
    }
    apply(params_mat, 1, sum_of_sq_dist)
}




# fit data and plot ------------------------------------------------------------

plot_probs <- function(prob_list) {

    plot_single <- function(df) {
        lines(x = df$x, y = df$y)
    }

    xvals <- unlist(purrr::map(prob_list, "x"))
    yvals <- unlist(purrr::map(prob_list, "y"))

    xlim <- c(min(xvals), max(xvals))
    ylim <- c(0, max(yvals))

    plot(-100, xlim = xlim, ylim = ylim, xlab = "Ovulation Day", ylab = "Standardized probability of conception")
    purrr::walk(prob_list, plot_single)
}




# toy examples -----------------------------------------------------------------

# # example skew normal distributions
# x <- seq(-5, 5, 0.01)
# plot(x, dskewnorm(x), type = "l", ylim = c(0, 0.7))
# points(x, dskewnorm(x, 0, 1, 1), type = "l", col = "red")
# points(x, dskewnorm(x, 0, 1, 2), type = "l", col = "blue")
# points(x, dskewnorm(x, 0, 1, 3), type = "l", col = "orange")
# points(x, dskewnorm(x, 0, 1, 4), type = "l", col = "green")

# # example skew normal distributions each normalized to have a maximum value of 1
# x <- seq(-5, 5, 0.01)
# plot(x, standardized_dskewnorm(x), type = "l", ylim = c(0, 1))
# points(x, standardized_dskewnorm(x, 0, 1, 1), type = "l", col = "red")
# points(x, standardized_dskewnorm(x, 0, 1, 2), type = "l", col = "blue")
# points(x, standardized_dskewnorm(x, 0, 1, 3), type = "l", col = "orange")
# points(x, standardized_dskewnorm(x, 0, 1, 4), type = "l", col = "green")
# points(x, standardized_dskewnorm(x, 0, 1, 5), type = "l", col = "purple")
# points(x, standardized_dskewnorm(x, 0, 1, 10), type = "l", col = "pink")




# find the value of nu ---------------------------------------------------------

find_probs <- function(decay) {
    0
}
