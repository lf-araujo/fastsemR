## helper-setup.R — run once before all tests

# Suppress OpenMx startup verbosity
if (requireNamespace("OpenMx", quietly = TRUE)) {
  suppressMessages({
    OpenMx::mxOption(key = "Number of Threads", value = 1L)
  })
  if (requireNamespace("umx", quietly = TRUE)) {
    suppressMessages(umx::umx_set_auto_run(FALSE))
  }
}
