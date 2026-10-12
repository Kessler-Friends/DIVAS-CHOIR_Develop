# Select a nuclear-norm penalty for an MCA matrix by holding out complete
# subject-by-variable responses (all indicator columns for one response).
categorical_nuclear_cv <- function(block, folds = 5L, lambda_grid = NULL,
                                   seed = NULL, max_iter = 200L, tol = 1e-5) {
  block <- as.data.frame(block)
  block[] <- lapply(block, factor)
  n <- nrow(block)
  if (ncol(block) < 2L || n < 3L) {
    stop("Nuclear CV needs at least two categorical variables and three subjects per block.")
  }
  if (length(folds) != 1L || !is.numeric(folds) || !is.finite(folds) ||
      folds != as.integer(folds) ||
      folds < 2L || folds > n) {
    stop("`categorical_cv_folds` must be an integer between 2 and the number of subjects.")
  }
  folds <- as.integer(folds)
  if (anyNA(block)) stop("Categorical blocks cannot contain missing responses.")

  # Preserve the caller's random-number state while making fold assignment
  # reproducible when a seed is supplied.
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  indicators <- lapply(block, function(x) {
    mm <- stats::model.matrix(~ x - 1)
    unname(mm)
  })
  widths <- vapply(indicators, ncol, integer(1))
  row_groups <- split(seq_len(sum(widths)), rep(seq_along(widths), widths))
  z <- t(do.call(cbind, indicators))
  full_p <- rowMeans(z)
  if (any(full_p == 0)) stop("Every categorical level must be observed at least once.")
  full_mca <- sweep(sweep(z, 1L, full_p, "-"), 1L, sqrt(full_p), "/")
  top <- svd(full_mca, nu = 0L, nv = 0L)$d[1L]
  if (!is.finite(top) || top <= 0) {
    stop("The categorical block has no MCA variation to estimate.")
  }
  if (is.null(lambda_grid)) lambda_grid <- top * exp(seq(log(0.02), log(0.8), length.out = 10L))
  if (!is.numeric(lambda_grid) || !length(lambda_grid) ||
      anyNA(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    stop("`categorical_lambda_grid` must contain positive finite numbers.")
  }
  lambda_grid <- sort(unique(lambda_grid))

  fold_id <- vapply(seq_along(indicators), function(j) {
    sample(rep(seq_len(folds), length.out = n))
  }, integer(n))
  if (is.null(dim(fold_id))) fold_id <- matrix(fold_id, ncol = 1L)
  errors <- matrix(NA_real_, nrow = folds, ncol = length(lambda_grid))

  for (f in seq_len(folds)) {
    observed <- matrix(FALSE, nrow(z), ncol(z))
    p_train <- numeric(nrow(z))
    for (j in seq_along(indicators)) {
      train <- fold_id[, j] != f
      rows <- row_groups[[j]]
      observed[rows, train] <- TRUE
      p_train[rows] <- rowSums(z[rows, train, drop = FALSE]) / sum(train)
    }
    present <- p_train > 0
    y <- sweep(z, 1L, p_train, "-")
    y[present, ] <- sweep(y[present, , drop = FALSE],
                          1L, sqrt(p_train[present]), "/")
    y[!present, ] <- 0
    for (i in seq_along(lambda_grid)) {
      h <- matrix(0, nrow(z), ncol(z))
      for (iter in seq_len(max_iter)) {
        next_input <- h
        next_input[observed] <- y[observed]
        s <- svd(next_input)
        d <- pmax(s$d - lambda_grid[i], 0)
        keep <- which(d > 0)
        h_new <- if (length(keep)) {
          (s$u[, keep, drop = FALSE] * rep(d[keep], each = nrow(z))) %*%
            t(s$v[, keep, drop = FALSE])
        } else matrix(0, nrow(z), ncol(z))
        if (sqrt(sum((h_new - h)^2)) <= tol * max(1, sqrt(sum(h^2)))) {
          h <- h_new
          break
        }
        h <- h_new
      }
      # A category absent from training has no defined MCA coordinate in
      # that fold, so it cannot contribute to held-out reconstruction error.
      held_out <- !observed & matrix(present, nrow(z), ncol(z))
      errors[f, i] <- mean((y[held_out] - h[held_out])^2)
    }
  }
  mean_error <- colMeans(errors)
  best <- which.min(mean_error)
  list(lambda = lambda_grid[best], cv = data.frame(
    lambda = lambda_grid, mean_error = mean_error,
    standard_error = apply(errors, 2L, stats::sd) / sqrt(folds)
  ))
}
