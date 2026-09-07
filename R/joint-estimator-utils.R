.joint_z_names <- function(data, z = NULL) {
  if (is.null(z)) {
    z <- grep("^Z[0-9]+$", names(data), value = TRUE)
  }
  if (!is.character(z) || length(z) == 0L || anyNA(z) || anyDuplicated(z)) {
    stop("`z` must identify at least one cheap covariate column.", call. = FALSE)
  }
  .check_required_columns(data, z)
  z
}

.joint_design <- function(data, first, z) {
  out <- as.matrix(data[, c(first, z), drop = FALSE])
  storage.mode(out) <- "double"
  out
}

.joint_formula <- function(response, terms, intercept = FALSE, environment = parent.frame()) {
  stats::as.formula(
    paste(response, "~", if (!intercept) "0 +" else "", paste(terms, collapse = " + ")),
    env = environment
  )
}

.joint_influence <- function(raw_score, information) {
  adjustment <- diag(1e-6, nrow(information))
  raw_score %*% t(solve(information + adjustment))
}

.joint_projection <- function(x, y) {
  cv <- glmnet::cv.glmnet(
    x = x, y = y, family = "mgaussian", intercept = FALSE, alpha = 0
  )
  model <- glmnet::glmnet(
    x = x, y = y, family = "mgaussian", intercept = FALSE,
    alpha = 0, lambda = cv$lambda.min
  )
  coef_matrix <- matrix(0, nrow = ncol(y), ncol = ncol(x))
  for (j in seq_len(ncol(y))) {
    coef_matrix[j, ] <- as.numeric(stats::coef(model)[[j]])[-1]
  }
  list(coef = coef_matrix, method = "glmnet")
}

.joint_solve_score <- function(initial, score) {
  as.numeric(nleqslv::nleqslv(initial, score)$x)
}

.joint_quadrature <- function(n_quad = 30L) {
  if (!is.numeric(n_quad) || length(n_quad) != 1L || is.na(n_quad) ||
      n_quad != as.integer(n_quad) || n_quad < 1L) {
    stop("`n_quad` must be a positive whole number.", call. = FALSE)
  }
  statmod::gauss.quad(as.integer(n_quad), kind = "hermite")
}

.joint_validate_nboot <- function(nboot) {
  if (!is.numeric(nboot) || length(nboot) != 1L || is.na(nboot) ||
      nboot != as.integer(nboot) || nboot < 0 || nboot == 1) {
    stop("`nboot` must be 0 or a whole number of at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

.joint_bootstrap_se <- function(data, nboot, fit_once, expected) {
  if (nboot == 0L) return(NULL)
  draws <- vector("list", nboot)
  completed <- 0L
  attempts <- 0L
  max_attempts <- nboot + max(10L, ceiling(nboot / 4))
  while (completed < nboot && attempts < max_attempts) {
    attempts <- attempts + 1L
    index <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
    draw <- try(fit_once(data[index, , drop = FALSE]), silent = TRUE)
    valid <- !inherits(draw, "try-error") &&
      identical(dim(draw), dim(expected)) && all(is.finite(draw))
    if (valid) {
      completed <- completed + 1L
      draws[[completed]] <- draw
    }
  }
  if (completed < 2L) {
    stop("Fewer than two joint bootstrap fits succeeded.", call. = FALSE)
  }
  if (completed < nboot) {
    warning(
      sprintf(
        "Only %d of %d requested joint bootstrap fits succeeded after %d attempts.",
        completed, nboot, attempts
      ),
      call. = FALSE
    )
  }
  draws <- draws[seq_len(completed)]
  array_draws <- array(
    unlist(draws, use.names = FALSE),
    dim = c(nrow(expected), ncol(expected), completed)
  )
  se <- apply(array_draws, c(1, 2), stats::sd)
  dimnames(se) <- dimnames(expected)
  attr(se, "bootstrap_attempts") <- attempts
  attr(se, "bootstrap_failures") <- attempts - completed
  se
}

.new_joint_fit <- function(call, family, estimates, se = NULL, details = list()) {
  if (is.null(se)) {
    se <- matrix(NA_real_, nrow(estimates), ncol(estimates), dimnames = dimnames(estimates))
  }
  results <- data.frame(
    method = rep(colnames(estimates), each = nrow(estimates)),
    term = rep(rownames(estimates), times = ncol(estimates)),
    estimate = as.vector(estimates),
    std.error = as.vector(se),
    row.names = NULL
  )
  structure(
    list(call = call, family = family, estimate = estimates, se = se,
         results = results, details = details),
    class = c("twostep2phase_joint", "twostep2phase")
  )
}

print.twostep2phase_joint <- function(x, ...) {
  cat("Joint two-step two-phase estimator\n")
  cat("Family:", x$family, "\n\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

.joint_fit_x_model <- function(data, x, x_star, z, weights,
                               type = c("linear", "gam")) {
  type <- match.arg(type)
  if (type == "linear") {
    formula <- .joint_formula(x, c(x_star, z), intercept = TRUE)
    fit <- stats::lm(formula, data = data, weights = weights)
  } else {
    smooth_terms <- paste0("s(", c(x_star, z), ")")
    formula_environment <- new.env(parent = asNamespace("gam"))
    formula_environment$weights <- weights
    formula <- stats::as.formula(
      paste(x, "~", paste(smooth_terms, collapse = " + ")),
      env = formula_environment
    )
    fit <- gam::gam(formula, family = stats::gaussian(), data = data,
                    weights = weights)
  }
  residuals <- if (type == "gam") fit$residuals else stats::residuals(fit)
  list(fit = fit, se = .weighted_resid_sd(residuals, weights))
}
