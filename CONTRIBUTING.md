# Contributing to fastsemR

Thanks for your interest in fastsemR.  This document describes how
to file issues, contribute code, and where to direct questions.

## Where to file an issue

There are two related repositories:

- **fastsemR** (this repo) — the R package, `umx` / OpenMx bridge,
  syntax translators, and R-side tests.  File issues here for:
  bugs in the R API, unexpected `summary()` output, problems with
  `run_fastsem()`, installation issues on Windows / macOS / Linux,
  documentation, vignettes.
  → <https://github.com/lf-araujo/fastsemR/issues>

- **fastsem** — the underlying Nim SEM engine.  File issues there
  for: numerical results, optimiser convergence, GPU / OpenCL
  kernels, performance, or the lavaan-style syntax parser.
  → <https://github.com/lf-araujo/fastsem/issues>

If you are not sure which side a bug lives on, file in fastsemR
and we will move it if needed.

## Reporting a bug

A useful bug report includes:

1. The output of `sessionInfo()` or `utils::sessionInfo()`.
2. The fastsem engine version reported by
   `getOption("fastsem.version")` after the package loads.
3. A minimal reproducible example — ideally one that fits in
   ~20 lines and uses simulated data or a built-in dataset.
4. The full error message and traceback (`traceback()` after the
   failure, or `options(error = recover)`).

## Submitting a pull request

1. Fork the repository and create a feature branch off `main`.
2. Run `R CMD check` locally; the GitHub Actions workflow runs the
   same checks on Linux, macOS, and Windows for every PR.
3. Add or update tests under `tests/testthat/` for any change in
   behaviour.  We use the testthat 3rd edition.
4. Update `NEWS.md` under "fastsemR (development version)" with a
   one-line user-facing summary of the change.
5. Open a PR against `main`.  Small, focused PRs are easier to
   review than large ones.

## Asking a question

For usage questions that are not bug reports, please open a GitHub
Discussion (or an issue with the `question` label) so the answer
is searchable for other users.

## Code of conduct

Please be respectful in all project communications.  Personal
attacks, harassment, and discriminatory language are not tolerated.
