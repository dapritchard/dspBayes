dskewnorm <- function(x, mu = 0, sigma = 1, shape = 0) {

  # standardize `x`
  x_star <- (x - mu) / sigma

  # calculate the density function
  2 / sigma * dnorm(x_star) * pnorm(shape * x_star)
}


standardized_dskewnorm <- function(x, mu = 0, sigma = 1, shape = 0) {
  norm_const <- find_skewnorn_mode(mu, sigma, shape)
  dskewnorm(x, mu, sigma, shape) / norm_const
}


find_skewnorn_mode <- function(mu = 0, sigma = 1, shape = 0) {
  f <- function(x) dskewnorm(x, mu, sigma, shape)
  if (sigma <= 0) return(0)
  optimize(f, c(-10, 10), maximum = TRUE)$objective
}


find_optimal_params <- function(x, y, starting = c(sigma = 1, shape = 0)) {

  f <- function(params){

    sigma = params[1L]
    shape = params[2L]

    y_hat <- standardized_dskewnorm(x, mu = 0, sigma = sigma, shape = shape)
    crossprod(y - y_hat)
  }

  optim(starting, f)
}


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




# historical data --------------------------------------------------------------

# prob_preg_list <- list(

#   columbo2006 = tibble::tibble(
#     x = c(   -8,    -7,    -6,    -5,    -4,    -3,    -2,    -1,     0,     1,     2,     3),
#     y = c(0.000, 0.012, 0.034, 0.035, 0.177, 0.249, 0.136, 0.228, 0.429, 0.114, 0.085, 0.033)
#   ),

#   stanford2007_19to26 = tibble::tibble(
#     x = c(    -8,     -7,     -6,     -5,     -4,     -3,     -2,     -1,      0,      1,      2),
#     y = c(0.0042, 0.0119, 0.0308, 0.0820, 0.2565, 0.2971, 0.5336, 0.3221, 0.1010, 0.0232, 0.0210)
#   ),

#   stanford2007_35to39 = tibble::tibble(
#     x = c(    -8,     -7,     -6,     -5,     -4,     -3,     -2,     -1,      0,      1,      2),
#     y = c(0.0039, 0.0108, 0.0265, 0.0658, 0.1728, 0.1951, 0.2901, 0.2029, 0.0798, 0.0205, 0.0187)
#   ),

#   lynch2006_bbt_colombo = tibble::tibble(
#     x = c(  -7,   -6,   -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3),
#     y = c(0.01, 0.03, 0.07, 0.18, 0.24, 0.26, 0.21, 0.10, 0.01, 0.04, 0.00)
#   ),

#   lynch2006_cm_columbo = tibble::tibble(
#     x = c(  -6,   -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3,    4),
#     y = c(0.00, 0.05, 0.08, 0.18, 0.11, 0.20, 0.18, 0.14, 0.07, 0.02, 0.01)
#   ),

#   lynch2006_lut_only_wilcox = tibble::tibble(
#     x = c(  -5,   -4,   -3,   -2,   -1,    0),
#     y = c(0.08, 0.17, 0.08, 0.36, 0.34, 0.36)
#   ),

#   lynch2006_lut_and_horm_dunson = tibble::tibble(
#     x = c(  -5,   -4,   -3,   -2,   -1,    0,    1),
#     y = c(0.00, 0.12, 0.12, 0.27, 0.24, 0.05, 0.00)
#   )
# )




# standardize data -------------------------------------------------------------

standardize_probs <- function(df) {
  max_idx <- which.max(df$y)
  df$x <- df$x - df$x[max_idx]
  df$y <- df$y / df$y[max_idx]
  df
}

# standardized_preg_list <- purrr::map(prob_preg_list, standardize_probs)
# standardized_preg_df <- do.call(rbind, standardized_preg_list)




# fit curve to historical data -------------------------------------------------

# fit <- grid_search_optimal_params(standardized_preg_df, seq(0.1, 5, 0.1), seq(-2, 2, 0.1)) %>%
#   dplyr::mutate(
#     sigma = purrr::map_dbl(param_vals, 1L),
#     shape = purrr::map_dbl(param_vals, 2L)
#   )

# conv_fit <- fit[fit$convergence, ] %>%
#   dplyr::select(grid_vals_sigma, grid_vals_shape, objective, sigma, shape)

# # TODO: think of a better way to pick a single value
# best_idx <- which(abs(conv_fit$objective - min(conv_fit$objective)) < 0.0001)
# best_idx <- best_idx[length(best_idx) %/% 2L]
# optimal_params <- c(
#   sigma = conv_fit$sigma[best_idx],
#   shape = conv_fit$shape[best_idx]
# )




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


# # non-standardized probabilities
# plot_probs(prob_preg_list)

# # standardized probabilities
# plot_probs(standardized_preg_list)
# x <- seq(-5, 5, 0.01)
# points(
#   x    = x,
#   y    = standardized_dskewnorm(x, 0, optimal_params["sigma"], optimal_params["shape"]),
#   type = "l",
#   col  = "red",
#   lwd  = 8
# )
