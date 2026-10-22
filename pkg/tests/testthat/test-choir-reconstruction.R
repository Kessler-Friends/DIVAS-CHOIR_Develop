test_that("fixed-score binary likelihood agrees with ordinary logistic regression", {
  x <- rep(c(-1, 1), each = 60)
  raw <- data.frame(a = factor(c(rep("a", 40), rep("b", 20), rep("a", 20), rep("b", 40))))
  V <- matrix(x + 3, ncol = 1L)
  fit <- choir_reconstruct(raw, V, "1", 1L)
  reference <- glm(I(raw$a == "b") ~ x, family = binomial())
  expect_equal(fit$nll, as.numeric(-logLik(reference)), tolerance = 1e-5)
  expect_equal(unname(fit$matLoadingMap[["1"]][2, 1] - fit$matLoadingMap[["1"]][1, 1]),
               unname(coef(reference)[2]), tolerance = 3e-4)
  expect_lt(fit$nll, fit$baseline_nll)
  expect_equal(rowMeans(fit$eta), c(0, 0), ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(fit$score_means, 3)
  expect_equal(fit$eta, fit$matBlockMap[["1"]])
  mca <- choir_reconstruct(MatCenterDIVASCHOIR(raw), V, "1", 1L)
  expect_equal(mca$probabilities, fit$probabilities, tolerance = 1e-6)
})

test_that("multiple categorical components add in logit space with variable-wise softmax", {
  x <- rep(c(-1, -1, 1, 1), each = 60)
  y <- rep(c(-1, 1, -1, 1), each = 60)
  a <- unlist(lapply(list(c(30,20,10), c(15,25,20), c(20,10,30), c(10,30,20)),
                     function(counts) rep(letters[1:3], counts)))
  raw <- data.frame(a = factor(a), b = factor(rep(c("x", "y"), 120)))
  fit <- choir_reconstruct(raw, cbind(x, y), c("3", "2"), c(1L, 1L))
  z <- MatCenterDIVASCHOIR(raw, method = "nuclear_norm")
  for (variable in attr(z, "divaschoir_variable_levels")) {
    rows <- match(variable$rows, rownames(z))
    expect_equal(sum(fit$mu[rows]), 0, tolerance = 1e-8)
    B <- do.call(cbind, fit$matLoadingMap)
    expect_equal(colSums(B[rows, , drop = FALSE]), c(0, 0), ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(colSums(fit$probabilities[rows, , drop = FALSE]), rep(1, 240),
                 ignore_attr = TRUE, tolerance = 1e-8)
  }
  expect_equal(Reduce(`+`, fit$matBlockMap), fit$eta)
  expect_equal(fit$theta, sweep(fit$eta, 1L, fit$mu, "+"))
  expect_equal(fit$residual, z - fit$probabilities, ignore_attr = TRUE)
  expect_equal(unname(rowMeans(fit$eta)), rep(0, nrow(z)), tolerance = 1e-8)
})

test_that("reconstruction routes categorical blocks to likelihood and continuous blocks to LS", {
  x <- rep(c(-1, 1), each = 60)
  raw <- data.frame(a = factor(c(rep("a",40),rep("b",20),rep("a",20),rep("b",40))))
  X <- rbind(x, 2*x + rep(c(-0.1,0.1),60))
  V <- matrix(x / sqrt(sum(x^2)), ncol = 1L)
  mca <- MatCenterDIVASCHOIR(raw)
  capture.output(fit <- DJIVEReconstructMJ(list(X, mca), c("continuous", "categorical"),
    list(`3`=V), list(`3`=1:2), "3", 0, categorical_blocks=list(NULL, raw)))
  basis <- fit$jointBasisMap[["3"]]
  expect_equal(as.numeric(fit$matLoadings[[1]][["3"]]), as.numeric(t(qr.solve(basis, t(X)))))
  catfit <- fit$categoricalReconstruction[[2]]
  expect_equal(fit$matBlocks[[2]][["3"]], catfit$eta)
  expect_true(all(is.finite(catfit$probabilities)))
  expect_null(fit$categoricalReconstruction[[1]])
  expect_equal(catfit$method, "multinomial_fixed_scores")
})

test_that("intercept-only categorical reconstruction returns marginal probabilities", {
  raw <- data.frame(a=factor(c("a","a","b","a","b","b")))
  fit <- choir_reconstruct(raw, matrix(0,6,0), character(), numeric())
  expect_equal(as.vector(fit$probabilities), rep(0.5,12))
  expect_length(fit$matBlockMap,0L)
  expect_equal(fit$eta, matrix(0,2,6,dimnames=dimnames(fit$eta)))
  expect_error(choir_reconstruct(raw, cbind(rep(1,6)), "1", 1L), "linearly independent")
})


test_that("main exposes likelihood coefficients and fixed scores for either extraction method", {
  x <- rep(c(-1, 1), each = 60)
  raw <- data.frame(a = factor(c(rep("a",40),rep("b",20),rep("a",20),rep("b",40))))
  V <- matrix(x/sqrt(sum(x^2)), ncol=1)
  local_mocked_bindings(
    DIVASCHOIRSignalExtract = function(datablock, method, ...) {
      prepared <- lapply(datablock, MatCenterDIVASCHOIR, method=method)
      list(preparedBlocks=prepared, analysisBlocks=prepared,
           VBars=list(V), UBars=list(matrix(c(-1,1)/sqrt(2),2,1)),
           phiBars=5, psiBars=5, EHats=list(matrix(0,2,120)), rBars=1,
           singVals=list(1), singValsHat=list(1), rSteps=list(1),
           VVHatCacheBars=list(list()), UUHatCacheBars=list(list()))
    },
    DJIVEJointStrucEstimateJP = function(...) list(outMap=list(`1`=V),
      keyIdxMap=list(`1`=1L), jointBlockOrder="1"),
    DJIVEAngleDiagnosticJP = function(...) NULL,
    .package="DIVAS")
  for (method in c("mca", "nuclear_norm")) {
    capture.output(result <- DIVASmain(categorical_datablock=list(cat=raw),
      categorical_method=method, nsim=2L, iprint=FALSE, ReturnDetail=TRUE))
    reconstruction <- result$categoricalReconstruction$cat
    expect_equal(unname(result$Loadings$cat), unname(reconstruction$matLoadingMap[["1"]]))
    expect_equal(unname(result$scoresList$cat), unname(reconstruction$basis))
    expect_lt(reconstruction$nll, reconstruction$baseline_nll)
  }
})
