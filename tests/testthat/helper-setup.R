## helper-setup.R — run once before all tests

# Suppress OpenMx startup verbosity
if (requireNamespace("OpenMx", quietly = TRUE)) {
  suppressMessages({
    OpenMx::mxOption(key = "Number of Threads", value = 1L)
  })
  if (requireNamespace("umx", quietly = TRUE)) {
    suppressMessages(umx::umx_set_auto_run(FALSE))
    # umx sets its variance-check defaults in .onAttach; under testthat umx
    # is only loadNamespace'd, so they must be set explicitly or umxACE()
    # crashes on getOption("umx_minVar") returning NULL.
    options(umx_minVar = 0.1, umx_maxVarRatio = 1000)
  }
}

# A local engine build (e.g. one carrying fixes not yet released) can be
# tested with FASTSEMR_LIB=/path/to/libfastsem_r.so.
if (nzchar(Sys.getenv("FASTSEMR_LIB")))
  suppressMessages(fastsemR::fastsem_load(Sys.getenv("FASTSEMR_LIB")))

skip_if_no_bridge_engine <- function() {
  skip_if_not(isTRUE(fastsemR:::.fastsem_env$loaded), "fastsem native library not loaded")
  ok <- tryCatch(fastsemR:::.require_engine_caps(fastsemR:::.bridge_caps),
                 error = function(e) FALSE)
  skip_if_not(isTRUE(ok), "loaded fastsem engine predates the run_fastsem() bridge")
}
