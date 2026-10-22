test_that("CHOIR preparation distinguishes observed MCA and latent centering", {
  raw <- data.frame(a = factor(c("a", "a", "a", "b")),
                    b = factor(c("x", "y", "x", "y")))
  mca <- MatCenterDIVASCHOIR(raw, method = "mca")
  nuclear <- MatCenterDIVASCHOIR(raw, method = "nuclear_norm")
  expect_true(all(nuclear %in% c(0, 1)))
  setup <- attr(nuclear, "divaschoir_setup")
  expect_equal(rowMeans(setup$eta), rep(0, nrow(nuclear)), ignore_attr = TRUE)
  for (rows in setup$groups) expect_equal(sum(setup$mu[rows]), 0, tolerance = 1e-12)
  expect_equal(unname(rowMeans(nuclear)), c(.75, .25, .5, .5))
  expect_equal(unname(rowMeans(mca)), rep(0, 4), tolerance = 1e-12)
  expected <- sweep(sweep(nuclear, 1L, rowMeans(nuclear), "-"),
                    1L, sqrt(rowMeans(nuclear)), "/")
  expect_equal(as.vector(mca), as.vector(expected))
  expect_error(MatCenterDIVASCHOIR(mca, method = "nuclear_norm"), "Prepared CHOIR block")
})

test_that("MCA estimation remains ordinary truncated SVD", {
  fit <- choir_estimate_mca(diag(c(12, 7, 3, 1)), 3L)
  expect_equal(fit$d, c(12, 7, 3))
  expect_equal(fit$signal, diag(c(12, 7, 3, 0)))
  expect_null(fit$lambda)
})

test_that("nuclear likelihood centers eta and identifies the category effects", {
  set.seed(51)
  raw <- as.data.frame(replicate(5, factor(sample(letters[1:3], 40, TRUE)),
                               simplify = FALSE))
  z <- MatCenterDIVASCHOIR(raw, method = "nuclear_norm")
  original <- z
  fit <- choir_estimate_nuclear_norm(z, lambda = 2)
  expect_true(fit$converged)
  expect_equal(rowMeans(fit$signal), rep(0, nrow(z)), ignore_attr = TRUE, tolerance = 1e-10)
  for (rows in attr(z, "divaschoir_setup")$groups) {
    expect_equal(sum(fit$mu[rows]), 0, tolerance = 1e-10)
    expect_equal(colSums(fit$signal[rows, ]), rep(0, ncol(z)),
                 ignore_attr = TRUE, tolerance = 1e-10)
    expect_equal(colSums(fit$probability_matrix[rows, ]), rep(1, ncol(z)), ignore_attr = TRUE)
  }
  expect_equal(z, original)
  expect_equal(as.vector(fit$theta), as.vector(sweep(fit$signal, 1L, fit$mu, "+")))
})

test_that("five-fold CV matches analytic held-out loss for intercept-only fits", {
  raw <- data.frame(a = factor(c("a", "a", "a", "a", "b", "b")),
                    b = factor(c("x", "x", "x", "y", "y", "y")))
  fit <- choir_estimate_nuclear_norm(raw, rank = 1L, lambda_grid = c(100, 200), tol = 1e-10)
  expected <- 0
  for (j in seq_along(raw)) for (fold in 1:5) {
    held <- fit$cv_folds[,j] == fold
    probabilities <- table(raw[[j]][!held]) / sum(!held)
    expected <- expected - sum(log(probabilities[as.character(raw[[j]][held])]))
  }
  expect_equal(fit$selection$cv_loss, rep(expected / (nrow(raw)*ncol(raw)), 2), tolerance = 1e-5)
  expect_equal(fit$rank, 0L)
  expect_identical(fit$selection_method, 'five_fold_response_cv')
  expect_equal(fit$lambda, fit$selection$lambda[which.min(fit$selection$cv_loss)])
})

test_that("response folds preserve training coverage and caller random state", {
  raw <- data.frame(a=factor(rep(c('a','b'),10)),
                    b=factor(rep(c('x','y'),10)),
                    c=factor(rep(c('l','m'),10)))
  set.seed(17)
  seed <- .Random.seed
  fit <- choir_estimate_nuclear_norm(raw,lambda_grid=100)
  expect_identical(.Random.seed,seed)
  expect_equal(sort(unique(as.vector(fit$cv_folds))),1:5)
  for(k in 1:5) {
    expect_true(all(rowSums(fit$cv_folds != k)>0))
    expect_true(all(colSums(fit$cv_folds != k)>0))
  }
  expect_equal(fit$selection$cv_loss,
               rowSums(fit$cv_fold_scores)/(nrow(raw)*ncol(raw)))
})

test_that("main prepares observations and passes the method to the CHOIR extractor", {
  raw <- data.frame(a = factor(rep(c("a", "b"), 10)),
                    b = factor(rep(c("x", "x", "y", "y"), 5)))
  for (choice in c("mca", "nuclear_norm")) {
    local_mocked_bindings(DIVASCHOIRSignalExtract = function(datablock, method, ...) {
      expect_identical(method, choice)
      expect_identical(datablock[[1L]], MatCenterDIVASCHOIR(raw, method=choice))
      stop("CHOIR route verified")
    }, .package = "DIVAS")
    expect_error(DIVASmain(categorical_datablock = list(cat = raw),
      categorical_method = choice, categorical_rank = 2), "CHOIR route verified")
  }
})

test_that("all nuclear refits receive unchanged one-hot observations", {
  set.seed(29)
  raw <- data.frame(a = factor(sample(c("a", "b"), 30, TRUE)),
                    b = factor(sample(c("x", "y"), 30, TRUE)),
                    c = factor(sample(c("l", "m"), 30, TRUE)))
  calls <- 0L
  penalties <- list()
  # Substitute only the numerical estimator; exercise real preparation,
  # probability generation, simulation and both bootstrap stages.
  local_mocked_bindings(choir_estimate_nuclear_norm = function(x, rank = NULL, lambda = NULL, lambda_grid = NULL, tol = 1e-7) {
    calls <<- calls + 1L
    penalties[calls] <<- list(lambda)
    if (!all(x %in% c(0, 1))) stop("Nuclear refit received transformed observations")
    mock_rank <- if (calls %in% c(3L, 1003L, 1004L, 1005L, 1006L)) 0L else 3L
    fit <- choir_estimate_mca(sweep(x, 1L, rowMeans(x), "-"), mock_rank)
    fit$method <- "nuclear_norm"
    fit$lambda <- if (is.null(lambda)) 0.7 else lambda
    fit$mu <- numeric(nrow(x))
    for (rows in attr(x, "divaschoir_setup")$groups) {
      a <- log(pmax(rowMeans(x[rows, , drop = FALSE]), 1e-8))
      fit$mu[rows] <- a - mean(a)
    }
    fit$probability_matrix <- x
    for (rows in attr(x, "divaschoir_setup")$groups) {
      logits <- sweep(fit$signal[rows, , drop = FALSE], 1L, fit$mu[rows], "+")
      a <- exp(sweep(logits, 2L, apply(logits, 2L, max), "-"))
      fit$probability_matrix[rows, ] <- sweep(a, 2L, colSums(a), "/")
    }
    fit
  }, .package = "DIVAS")
  capture.output(fit <- DIVASCHOIRSignalExtract(list(cat = raw), fixed_rank = 2L,
    nsim = 2L, method = "nuclear_norm", seed = 29))
  expect_equal(calls, 1006L)
  expect_equal(fit$estimation[[1L]]$rank, 2L)
  expect_equal(fit$estimation[[1L]]$untruncated_rank, 3L)
  expect_equal(fit$estimation[[1L]]$null_reference_rank, 2L)
  expect_true(all(fit$estimation[[1L]]$null_ranks %in% c(0L, 1L)))
  expect_null(penalties[[1L]])
  expect_true(all(vapply(penalties[-1L], identical, logical(1), 0.7)))
  expect_true(all(fit$preparedBlocks[[1L]] %in% c(0, 1)))
  expect_equal(rowMeans(fit$analysisBlocks[[1L]]), rep(0, nrow(fit$analysisBlocks[[1L]])),
               ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(fit$estimation[[1L]]$angle_ranks, c(0L, 0L))
  expect_equal(fit$estimation[[1L]]$retained_ranks, c(0L, 0L))
  expect_true(is.na(fit$nullAngleCaches[[1L]]$subject[1L]))
  expect_equal(fit$phiBars, 90)
  expect_equal(fit$VVHatCacheBars[[1L]][[2L]], matrix(0, 2, 2))
})


test_that("CV-selected zero rank returns empty spaces without forced directions", {
  raw <- data.frame(a = factor(rep(c("a", "b"), each = 3)),
                    b = factor(rep(c("x", "y"), 3)))
  capture.output(fit <- DIVASCHOIRSignalExtract(list(cat = raw), nsim = 2L,
    method = "nuclear_norm", lambda_grid = 100))
  expect_equal(fit$rBars, 0)
  expect_equal(ncol(fit$VBars[[1L]]), 0L)
  expect_equal(ncol(fit$UBars[[1L]]), 0L)
  expect_identical(fit$estimation[[1L]]$selection_method, 'five_fold_response_cv')
  expect_length(fit$estimation[[1L]]$angle_ranks, 0L)
})

test_that("nuclear culling can reject every direction and skips retained refits", {
  set.seed(19)
  raw <- data.frame(a = factor(sample(c("a", "b"), 30, TRUE)),
                    b = factor(sample(c("x", "y"), 30, TRUE)),
                    c = factor(sample(c("l", "m"), 30, TRUE)))
  penalties <- list()
  # Exercise real null simulation, angle bootstrapping, culling and empty-space
  # outputs; replace only the expensive nuclear optimization.
  local_mocked_bindings(choir_estimate_nuclear_norm = function(x, rank = NULL,
      lambda = NULL, lambda_grid = NULL, tol = 1e-7) {
    penalties[length(penalties) + 1L] <<- list(lambda)
    fit <- choir_estimate_mca(sweep(x, 1L, rowMeans(x), "-"), 3L)
    fit$lambda <- if (is.null(lambda)) 0.7 else lambda
    fit$mu <- numeric(nrow(x))
    fit$probability_matrix <- x
    fit
  }, .package = "DIVAS")
  capture.output(fit <- DIVASCHOIRSignalExtract(list(cat = raw), fixed_rank = 2L,
    nsim = 2L, method = "nuclear_norm", cull_rank = TRUE, cull = 0, seed = 19))
  expect_length(penalties, 1004L) # Original, null reference, 1,000 null, 2 angle fits.
  expect_null(penalties[[1L]])
  expect_true(all(vapply(penalties[-1L], identical, logical(1), 0.7)))
  expect_equal(fit$estimation[[1L]]$rank, 2L)
  expect_equal(fit$rBars, 0)
  expect_false(any(fit$rankDiagnostics[[1L]]$retained))
  expect_equal(dim(fit$VBars[[1L]]), c(nrow(raw), 0L))
  expect_equal(dim(fit$UBars[[1L]]), c(6L, 0L))
  expect_equal(fit$phiBars, 90)
  expect_equal(fit$psiBars, 90)
  expect_length(fit$VVHatCacheBars[[1L]], 0L)
  expect_length(fit$UUHatCacheBars[[1L]], 0L)
  expect_length(fit$estimation[[1L]]$retained_ranks, 0L)
  expect_length(fit$estimation[[1L]]$retained_lambdas, 0L)
  expect_length(fit$estimation[[1L]]$angle_ranks, 2L)
  expect_length(fit$estimation[[1L]]$null_ranks, 1000L)
})


test_that("five-fold CV selects a nonzero signal for strongly related binary traits", {
  raw <- data.frame(a = factor(rep(c("a", "b"), each = 4)),
                    b = factor(rep(c("x", "y"), each = 4)),
                    c = factor(rep(c("l", "m"), each = 4)))
  fit <- choir_estimate_nuclear_norm(raw, rank = 3L, lambda_grid = c(0.3, 1, 3),
                                    max_iter = 2000L)
  expect_true(all(is.finite(fit$selection$cv_loss)))
  expect_equal(fit$lambda, fit$selection$lambda[which.min(fit$selection$cv_loss)])
  expect_equal(fit$rank, 1L)
  expect_equal(fit$objective, -sum(fit$log_probs[MatCenterDIVASCHOIR(raw,
    method = "nuclear_norm") == 1]) + fit$lambda * sum(fit$d), tolerance = 1e-6)
  expect_error(suppressWarnings(choir_estimate_nuclear_norm(raw,
                 lambda_grid = c(0.3, 1), max_iter = 1L)),
               "Nuclear fit did not converge")
})

test_that("main accepts a CV-selected intercept-only categorical block", {
  raw <- data.frame(a = factor(rep(c("a", "b"), each = 3)),
                    b = factor(rep(c("x", "y"), 3)))
  expect_warning(capture.output(result <- DIVASmain(categorical_datablock = list(cat = raw),
    categorical_method = "nuclear_norm", categorical_lambda_grid = 100,
    nsim = 2L, iprint = FALSE, ReturnDetail = TRUE)), NA)
  expect_equal(result$Phase1_categorical$rBars, 0)
  expect_null(result$diagnostic_plots)
})


test_that("unculled nuclear runs can skip null fits and use a separate refit tolerance", {
  set.seed(19)
  raw <- data.frame(a=factor(sample(c("a","b"),30,TRUE)),
                    b=factor(sample(c("x","y"),30,TRUE)),
                    c=factor(sample(c("l","m"),30,TRUE)))
  tolerances <- numeric()
  local_mocked_bindings(choir_estimate_nuclear_norm=function(x, rank=NULL,
      lambda=NULL, lambda_grid=NULL, tol=1e-7) {
    tolerances <<- c(tolerances,tol)
    fit <- choir_estimate_mca(sweep(x,1L,rowMeans(x),"-"),3L)
    fit$lambda <- 0.7; fit$mu <- numeric(nrow(x)); fit$probability_matrix <- x
    fit
  }, .package="DIVAS")
  output <- capture.output(fit <- DIVASCHOIRSignalExtract(list(cat=raw),
    fixed_rank=2L,nsim=2L,method="nuclear_norm",simulate_null=FALSE,
    bootstrap_tol=1e-5,seed=19))
  expect_equal(tolerances,c(1e-7,rep(1e-5,4)))
  expect_length(fit$estimation[[1]]$null_ranks,0L)
  expect_equal(fit$estimation[[1]]$angle_ranks,c(2L,2L))
  expect_equal(fit$estimation[[1]]$retained_ranks,c(2L,2L))
  expect_true(any(grepl("Skipping categorical null",output)))
  expect_error(DIVASCHOIRSignalExtract(list(cat=raw),cull_rank=TRUE,
    simulate_null=FALSE),"requires null")
})
