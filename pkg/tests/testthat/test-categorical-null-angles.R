test_that("loading diagnostics use retained-space angles and omit degenerate null bounds", {
  caches <- lapply(c(10, 20, 30), function(angle) diag(c(1, cos(angle * pi / 180))))
  bounds <- choir_loading_diagnostic_bounds(caches, 0)
  expect_equal(bounds$perturbation, unname(quantile(c(10, 20, 30), 0.95)))
  expect_true(is.na(bounds$random))
  expect_equal(choir_loading_diagnostic_bounds(caches, 75)$random, 75)
  expect_equal(choir_loading_diagnostic_bounds(list(matrix(0, 2, 2)), 75)$perturbation, 90)
  expect_equal(choir_loading_diagnostic_bounds(list(diag(1 + 1e-12, 2)), 75)$perturbation, 0)
})

test_that("categorical null directions are compared to one fixed null subspace", {
  raw <- data.frame(
    a = factor(rep(c("a", "a", "a", "b", "c", "c"), 10)),
    b = factor(rep(c("x", "y", "y", "y", "x"), 12)),
    c = factor(rep(c("l", "m", "n", "l"), 15))
  )
  block <- MatCenterDIVASCHOIR(raw)
  metadata <- attr(block, "divaschoir_variable_levels")
  p_hat <- attr(block, "divaschoir_p_hat")

  # Reproduce multiple directions against the same reference subspace.
  draw_null <- function() {
    draw <- lapply(metadata, function(variable) {
      p <- p_hat[match(variable$rows, rownames(block))]
      mu <- log(p) - mean(log(p))
      probabilities <- exp(mu - max(mu))
      probabilities <- probabilities / sum(probabilities)
      expect_equal(probabilities, p)
      factor(vapply(seq_len(ncol(block)), function(i) {
        sample(variable$levels, 1L, prob = probabilities)
      }, character(1)), levels = variable$levels)
    })
    names(draw) <- vapply(metadata, `[[`, character(1), "variable")
    svd(MatCenterDIVASCHOIR(as.data.frame(draw), allow_empty_levels = TRUE),
        nu = 2L, nv = 2L)
  }
  angle <- function(direction, basis) {
    projection <- basis %*% crossprod(basis, direction)
    cosine <- sqrt(sum(projection^2)) / sqrt(sum(direction^2))
    acos(pmin(1, pmax(0, cosine))) * 180 / pi
  }
  set.seed(123)
  reference <- draw_null()
  expected_subject <- expected_feature <- numeric(3L)
  for (i in seq_len(3L)) {
    direction <- draw_null()
    expected_subject[i] <- angle(direction$v[, 1L], reference$v)
    expected_feature[i] <- angle(direction$u[, 1L], reference$u)
  }

  fit <- NULL
  capture.output(fit <- DIVASCHOIRSignalExtract(
    list(categorical = block), fixed_rank = 2L, nsim = 4L, seed = 123
  ))
  cache <- fit$nullAngleCaches[[1L]]
  expect_length(cache$subject, 1000L)
  expect_length(cache$feature, 1000L)
  expect_equal(cache$subject[1:3], expected_subject)
  expect_equal(cache$feature[1:3], expected_feature)
  expect_true(all(is.finite(unlist(cache))))
  expect_true(all(unlist(cache) >= 0 & unlist(cache) <= 90))
  diagnostics <- fit$rankDiagnostics[[1L]]
  expect_equal(diagnostics$subject_random_deg,
               rep(unname(quantile(cache$subject, 0.05)), 2L))
  expect_equal(diagnostics$feature_random_deg,
               rep(unname(quantile(cache$feature, 0.05)), 2L))
  expect_equal(diagnostics$subject_cutoff_deg,
               diagnostics$subject_random_deg * 0.5)
  expect_equal(fit$rBars, 2)

  culled <- NULL
  capture.output(culled <- DIVASCHOIRSignalExtract(
    list(categorical = block), fixed_rank = 2L, nsim = 4L, seed = 123,
    cull_rank = TRUE
  ))
  expect_equal(culled$nullAngleCaches, fit$nullAngleCaches)
  expected_retained <- diagnostics$subject_q95_deg < diagnostics$subject_cutoff_deg
  expect_equal(culled$rankDiagnostics[[1L]]$retained, expected_retained)
  expect_equal(culled$rBars, sum(expected_retained))
})

test_that("main accepts a categorical block culled to zero rank", {
  raw <- data.frame(
    a = factor(rep(c("a", "a", "a", "b", "c", "c"), 10)),
    b = factor(rep(c("x", "y", "y", "y", "x"), 12)),
    c = factor(rep(c("l", "m", "n", "l"), 15))
  )
  expect_warning(capture.output(fit <- DIVASmain(
    categorical_datablock = list(cat = raw), categorical_method = "mca",
    categorical_rank = 2L, categorical_cull = TRUE, nsim = 4L, seed = 123,
    iprint = FALSE, ReturnDetail = TRUE
  )), NA)
  expect_equal(fit$Phase1_categorical$estimation[[1L]]$rank, 2L)
  expect_equal(fit$rBars, 0)
  expect_equal(ncol(fit$VBars[[1L]]), 0L)
  expect_equal(ncol(fit$UBars[[1L]]), 0L)
  expect_null(fit$Scores)
  expect_null(fit$diagnostic_plots)
})

test_that("main forwards the null count separately from bootstrap count", {
  raw <- data.frame(a = factor(rep(c("a", "b"), 20)),
                    b = factor(rep(c("x", "x", "y", "y"), 10)))
  output <- capture.output(fit <- DIVASmain(
    categorical_datablock = list(cat = raw), categorical_method = "mca",
    categorical_rank = 2L, categorical_cull = TRUE, categorical_null_nsim = 7L,
    categorical_cull_fraction = 0.5,
    nsim = 3L, seed = 123, iprint = FALSE, ReturnDetail = TRUE
  ))
  categorical <- fit$Phase1_categorical
  expect_length(categorical$nullAngleCaches[[1L]]$subject, 7L)
  expect_length(categorical$nullAngleCaches[[1L]]$feature, 7L)
  expect_length(categorical$estimation[[1L]]$null_ranks, 7L)
  expect_equal(categorical$estimation[[1L]]$null_nsim, 7L)
  expect_length(categorical$estimation[[1L]]$angle_ranks, 3L)
  expect_equal(categorical$rankDiagnostics[[1L]]$subject_cutoff_deg,
               categorical$rankDiagnostics[[1L]]$subject_random_deg * 0.5)
  expect_true(any(grepl("Null simulation: 7/7", output, fixed = TRUE)))
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "200")) {
    expect_error(DIVASCHOIRSignalExtract(list(cat = raw), null_nsim = bad),
                 "null_nsim must be a positive integer")
  }
  for (bad in list(0, -1, 0.51, NA_real_, Inf, c(0.4, 0.5), "0.5")) {
    expect_error(DIVASmain(categorical_cull_fraction = bad),
                 "categorical_cull_fraction must be a number", fixed = TRUE)
  }
})
