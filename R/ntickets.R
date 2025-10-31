#'Optimal tickets sales
#'
#' @param N The number of seats available
#' @param gamma Probability of overbooking
#' @param p Probability a passenger with a ticket shows
#'
#' @returns Two plots and a list containing five components
#' @export
#'
#' @examples
#' ntickets(N=400,gamma = 0.02, p = 0.95)
ntickets <- function(N, gamma, p) {
  # Candidate ticket values
  nvals <- seq(N, N + 50)

  # Objective functions
  obj_discrete <- function(n) abs(1 - gamma - stats::pbinom(N, size = n, prob = p))
  obj_normal <- function(n) {
    mu <- n * p
    sigma <- sqrt(n * p * (1 - p))
    abs(1 - gamma - stats::pnorm(N + 0.5, mean = mu, sd = sigma))
}

  # Evaluate over range
  f_discrete <- sapply(nvals, obj_discrete)
  f_normal <- sapply(nvals, obj_normal)

  # ticket counts
  nd <- nvals[which.min(f_discrete)]
  nc <- nvals[which.min(f_normal)]

  par(mfrow = c(1, 2))  # one window, two plots side by side

  plot(nvals, f_discrete, type = "l", col = "blue", lwd = 2,
       ylab = "Objective Value", xlab = "n",
       main = "Discrete Objective")
  graphics::abline(h = 0, col = "gray", lty = 3)

  plot(nvals, f_normal, type = "l", col = "red", lwd = 2,
       ylab = "Objective Value", xlab = "n",
       main = "Normal Approximation")
  graphics::abline(h = 0, col = "gray", lty = 3)


  # Return named list
  return(list(nd = nd, nc = nc, N = N, p = p, gamma = gamma))
}

