## install.R — public binary management functions

#' Install or re-install the fastsem native binary.
#'
#' Downloads the pre-built shared library for your platform from the fastsem
#' GitHub Releases page and loads it into the current R session.  The binary
#' is cached in `tools::R_user_dir("fastsemR", "cache")` and re-used on
#' subsequent calls (including package startup via [`.onLoad`]).
#'
#' Platform detection is automatic:
#' * Linux x86-64  → `libfastsem_r_linux_x86_64.so`
#' * macOS x86-64  → `libfastsem_r_macos_x86_64.dylib`
#' * macOS arm64   → `libfastsem_r_macos_arm64.dylib`
#'
#' @param force Logical (default `FALSE`).  When `TRUE`, re-download even if a
#'   cached binary already exists.  Use after a new fastsem release.
#'
#' @return Invisibly, the absolute path to the loaded binary.
#'
#' @seealso [fastsem_update()] for a one-liner force-download;
#'   [fastsem_load()] to use a locally compiled binary.
#'
#' @export
#' @examples
#' \dontrun{
#' # First-time install (or after a fresh R install)
#' fastsem_install()
#'
#' # Upgrade to the latest release
#' fastsem_install(force = TRUE)
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
#' A convenience shorthand for `fastsem_install(force = TRUE)`.  Use this
#' whenever a new version of fastsem is released to pull the updated binary.
#'
#' @return Invisibly, the absolute path to the loaded binary.
#'
#' @seealso [fastsem_install()] for the full install function.
#'
#' @export
#' @examples
#' \dontrun{
#' fastsem_update()
#' }
fastsem_update <- function() fastsem_install(force = TRUE)

#' Load a locally built fastsem shared library.
#'
#' Use this when you have compiled the fastsem library yourself (e.g. with
#' `nimble buildRnim` in the fastsem source tree) and want to use that build
#' instead of the downloaded release binary.  Also useful for development and
#' testing without internet access.
#'
#' @param lib_path Character.  Path to:
#'   * `libfastsem_r.so` on Linux,
#'   * `libfastsem_r.dylib` on macOS, or
#'   * `libfastsem_r.dll` on Windows.
#'
#'   When `NULL` (default) the function looks for any of those file names in
#'   the current working directory.
#'
#' @return Invisibly `TRUE` on success.
#'
#' @seealso [fastsem_install()] to download a pre-built release binary.
#'
#' @export
#' @examples
#' \dontrun{
#' # Point to a local build
#' fastsem_load("~/fastsem/libfastsem_r.so")
#'
#' # Load from the current working directory (default search)
#' fastsem_load()
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
