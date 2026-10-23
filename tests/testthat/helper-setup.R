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
