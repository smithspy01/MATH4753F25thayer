#' Birthday
#'
#' @param x vector
#'
#' @returns a function of x
#' @export
#'
#' @examples
#' birthday(20:24)
birthday <- function(x){
  1 - exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}
