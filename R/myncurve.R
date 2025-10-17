#' List
#'
#' @param mu mean
#' @param sigma standard deviation
#' @param a area under the curve of which is shaded up to
#'
#' @returns a list with mu, sigma, and probability of a
#' @export
#'
#' @examples
#' myncurve(mu = 10, sigma = 5, a = 6)

myncurve <- function(mu, sigma, a) {
  curve(dnorm(x, mean = mu, sd = sigma),
        xlim = c(mu - 3 * sigma, mu + 3 * sigma),
        main = paste("P(X <=", a, ")"),
        ylab = "Density")

  x_vals <- seq(mu - 3 * sigma, a, length = 1000)
  y_vals <- dnorm(x_vals, mean = mu, sd = sigma)

  polygon(c(x_vals, a), c(y_vals, 0), col = "skyblue")

  prob <- pnorm(a, mean = mu, sd = sigma)

  list(mu = mu, sigma = sigma, prob = prob)
}
