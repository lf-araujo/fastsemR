## install.R — public binary management functions

#' Install or re-install the fastsem native binary.
#'
#' Downloads the pre-built shared library for your platform from the fastsem
#' GitHub Releases page and loads it into the current R session.  By default
#' the binary is cached in `tools::R_user_dir("fastsemR", "cache")` and
#' re-used on subsequent calls (including package startup).
#'
#' @param force Logical.  If `TRUE`, re-download even if a cached binary
#'   already exists.  Useful after a new fastsem release.
#'
#' @return Invisibly, the path to the loaded binary.
#' @export
#' @examples
#' \dontrun{
#' fastsem_install()          # first-time install
#' fastsem_install(force = TRUE)  # upgrade to latest release
#' }
fastsem_install <- function(force = FALSE) {
  if (isTRUE(.fastsem_env$loaded) && !force) {
    message("fastsem: already loaded from ", .fastsem_env$lib_path,
            "\n  Use fastsem_install(force = TRUE) to re-download.")
    return(invisible(.fastsem_env$lib_path))
  }

  lib_path <- .fastsem_download(force = force)

  if (isTRUE(.fastsem_env$loaded) && !is.null(.fastsem_env$lib_path)) {
    tryCatch(dyn.unload(.fastsem_env$lib_path), error = function(e) NULL)
    .fastsem_env$loaded     <- FALSE
    .fastsem_env$lib_path   <- NULL
    .fastsem_env$dll_handle <- NULL
  }

  handle <- dyn.load(lib_path)
  .fastsem_env$loaded     <- TRUE
  .fastsem_env$lib_path   <- lib_path
  .fastsem_env$dll_handle <- handle
  message("fastsem: loaded from ", lib_path)
  invisible(lib_path)
}

#' Force-download the latest fastsem binary and reload.
#'
#' Shorthand for `fastsem_install(force = TRUE)`.
#'
#' @return Invisibly, the path to the loaded binary.
#' @export
#' @examples
#' \dontrun{
#' fastsem_update()
#' }
fastsem_update <- function() fastsem_install(force = TRUE)

#' Load a locally built fastsem shared library.
#'
#' Use this when you have compiled the library yourself (e.g. with
#' `nimble buildRnim` in the fastsem source tree) and want to use that
#' build instead of the downloaded release binary.
#'
#' @param lib_path Path to `libfastsem_r.so` (Linux),
#'   `libfastsem_r.dylib` (macOS), or `libfastsem_r.dll` (Windows).
#'   When `NULL` (default) the function looks in the current working
#'   directory.
#'
#' @return Invisibly `TRUE` on success.
#' @export
#' @examples
#' \dontrun{
#' fastsem_load("~/fastsem/libfastsem_r.so")
#' }
fastsem_load <- function(lib_path = NULL) {
  if (is.null(lib_path)) {
    candidates <- c(
      file.path(getwd(), "libfastsem_r.so"),
      file.path(getwd(), "libfastsem_r.dylib"),
      file.path(getwd(), "libfastsem_r.dll")
    )
    lib_path <- Filter(file.exists, candidates)[1]
    if (length(lib_path) == 0 || is.na(lib_path))
      stop("fastsem: cannot find libfastsem_r.{so,dylib,dll} in ",
           getwd(), ".\n  Pass the path explicitly or run fastsem_install().")
  }

  if (isTRUE(.fastsem_env$loaded) && !is.null(.fastsem_env$lib_path)) {
    if (normalizePath(.fastsem_env$lib_path, mustWork = FALSE) ==
        normalizePath(lib_path, mustWork = FALSE)) {
      message("fastsem: already loaded from ", .fastsem_env$lib_path)
      return(invisible(TRUE))
    }
    tryCatch(dyn.unload(.fastsem_env$lib_path), error = function(e) NULL)
    .fastsem_env$loaded     <- FALSE
    .fastsem_env$dll_handle <- NULL
  }

  handle <- dyn.load(lib_path)
  .fastsem_env$loaded     <- TRUE
  .fastsem_env$lib_path   <- lib_path
  .fastsem_env$dll_handle <- handle
  message("fastsem: loaded from ", lib_path)
  invisible(TRUE)
}
