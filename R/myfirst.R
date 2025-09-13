#' My first function
#'
#' @param x A numeric vector.
#'
#' @returns A list with components `x` and `y`, where `y` is the square of `x`
#' @export
#'
#' @examples
#' myfirstfunction(1:10)
myfirstfunction <- function(x) {
  y <- x^2
  Plot(y~x)
  list(x = x, y = y)
}
