##' Print a COMPASSResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{COMPASS}} call.
##'
##' @param x An object of class \code{COMPASSResult}.
##' @param ... Optional arguments; currently unused.
##' @export
##' @examples
##' print(CR)
print.COMPASSResult <- function(x, ...) {
  n <- nrow(x$data$n_s)
  cat("A COMPASS model fit on", n, "paired samples.\n")
}

##' Summarize a COMPASSResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{COMPASS}} call.
##'
##' @param object An object of class \code{COMPASSResult}.
##' @param ... Optional arguments; currently unused.
##' @export
##' @examples
##' print(CR)
summary.COMPASSResult <- function(object, ...) {
  n <- nrow(object$data$n_s)
  cat("A COMPASSResult model fit on", n, "paired samples.\n")
}
