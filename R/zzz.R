## zzz.R — package hooks and binary lifecycle management

.fastsem_env <- new.env(parent = emptyenv())
.fastsem_env$loaded     <- FALSE
.fastsem_env$lib_path   <- NULL
.fastsem_env$dll_handle <- NULL   # DLLInfo returned by dyn.load()

# ── Platform detection ─────────────────────────────────────────────────────────

.fastsem_lib_filename <- function() {
  os   <- Sys.info()[["sysname"]]
  arch <- Sys.info()[["machine"]]   # "x86_64" or "arm64" / "aarch64"
  arch_tag <- if (grepl("arm|aarch", arch, ignore.case = TRUE)) "arm64" else "x86_64"
  win_arch  <- if (grepl("arm|aarch", arch, ignore.case = TRUE)) "arm64" else "x64"
  switch(os,
    Linux   = sprintf("libfastsem_r-linux-%s.so",    arch_tag),
    Darwin  = sprintf("libfastsem_r-macos-%s.dylib", arch_tag),
    Windows = sprintf("libfastsem_r-windows-%s.dll", win_arch),
    stop("fastsem: unsupported OS '", os, "'")
  )
}

# ── Cache directory ────────────────────────────────────────────────────────────

.fastsem_cache_dir <- function() {
  d <- tryCatch(
    tools::R_user_dir("fastsemR", "cache"),
    error = function(e) file.path(Sys.getenv("HOME"), ".cache", "fastsemR")
  )
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# ── Download binary from GitHub Releases ──────────────────────────────────────

.fastsem_download <- function(force = FALSE) {
  filename  <- .fastsem_lib_filename()
  cache_dir <- .fastsem_cache_dir()
  dest      <- file.path(cache_dir, filename)

  if (!force && file.exists(dest)) return(dest)

  base_url <- "https://github.com/lf-araujo/fastsemR/releases/download/0.1"
  url      <- paste0(base_url, "/", filename)

  message("fastsem: downloading binary\n  ", url)
  utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  dest
}

# ── .onLoad ────────────────────────────────────────────────────────────────────

.onLoad <- function(libname, pkgname) {
  tryCatch({
    lib_path <- .fastsem_download(force = FALSE)
    handle   <- dyn.load(lib_path)
    .fastsem_env$loaded     <- TRUE
    .fastsem_env$lib_path   <- lib_path
    .fastsem_env$dll_handle <- handle
  }, error = function(e) {
    packageStartupMessage(
      "fastsem: native library not loaded.\n",
      "  Run fastsem_install() to download the pre-built binary, or\n",
      "  fastsem_load(\"/path/to/libfastsem_r.so\") for a local build.\n",
      "  Details: ", conditionMessage(e)
    )
  })
}

.onUnload <- function(libpath) {
  if (isTRUE(.fastsem_env$loaded) && !is.null(.fastsem_env$lib_path)) {
    tryCatch(dyn.unload(.fastsem_env$lib_path), error = function(e) NULL)
    .fastsem_env$loaded     <- FALSE
    .fastsem_env$lib_path   <- NULL
    .fastsem_env$dll_handle <- NULL
  }
}
