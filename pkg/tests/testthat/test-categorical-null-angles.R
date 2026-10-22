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
               diagnostics$subject_random_deg * 0.382)
  expect_equal(fit$rBars, 2)

  culled <- NULL
  capture.output(culled <- DIVASCHOIRSignalExtract(
    list(categorical = block), fixed_rank = 2L, nsim = 4L, seed = 123,
    cull_rank = TRUE
  ))
  expect_equal(culled$nullAngleCaches, fit$nullAngleCaches)
  expected_retained <- diagnostics$subject_q95_deg < diagnostics$subject_cutoff_deg
  expected_retained[1L] <- TRUE
  expect_equal(culled$rankDiagnostics[[1L]]$retained, expected_retained)
  expect_equal(culled$rBars, sum(expected_retained))
})
