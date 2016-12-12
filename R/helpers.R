#' Test if a package is installed.
#'
#'
#' @param pkg Package name to test.
TestRequiredPkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(pkg, "needed for this function to work. Please install it."),
         call. = FALSE)
  }
}

#' Print a message if getOption("tess3r.debug") is not null
#'
#'
#' @param msg The message
DebugMessage <- function(msg) {
  if (!is.null(getOption("tess3r.debug"))) {
    message("= ", msg)
  }
}
