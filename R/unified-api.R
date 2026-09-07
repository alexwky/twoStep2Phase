.analysis_z_names <- function(data, z) {
  if (is.null(z)) return(character())
  if (!is.character(z) || anyNA(z) || anyDuplicated(z)) {
    stop("`z` must be a character vector of unique column names.", call. = FALSE)
  }
  .check_required_columns(data, z)
  z
}

.analysis_stratum <- function(data, stratum) {
  if (!is.null(stratum)) {
    if (!is.character(stratum) || length(stratum) != 1L || is.na(stratum)) {
      stop("`stratum` must be NULL or one column name.", call. = FALSE)
    }
    .check_required_columns(data, stratum)
    return(list(data = data, stratum = stratum, default = FALSE))
  }
  name <- ".twostep2phase_stratum"
  while (name %in% names(data)) name <- paste0(name, "_")
  data[[name]] <- 1L
  list(data = data, stratum = name, default = TRUE)
}

.analysis_x_models <- function(z, x_models) {
  if (is.null(x_models) || length(x_models) == 0L) return(character())
  if (!is.character(x_models) || anyNA(x_models)) {
    stop("`x_models` must be NULL or a character vector.", call. = FALSE)
  }
  x_models <- unique(x_models)
  if (length(z) == 0L) {
    if (!identical(x_models, "kernel")) {
      stop("Without `Z`, `x_models` must be \"kernel\" or NULL.", call. = FALSE)
    }
    return(x_models)
  }
  if (any(!x_models %in% c("linear", "gam"))) {
    stop(
      "With `Z`, `x_models` may contain only \"linear\" and/or \"gam\", or be NULL.",
      call. = FALSE
    )
  }
  .joint_x_models(x_models)
}

.analysis_initial <- function(initial, z) {
  if (is.null(initial)) return(NULL)
  if (!is.numeric(initial) || anyNA(initial) || any(!is.finite(initial))) {
    stop("`initial` must be NULL or a finite numeric vector.", call. = FALSE)
  }
  if (length(z) == 0L) {
    if (length(initial) != 2L || initial[1L] >= initial[2L]) {
      stop(
        "Without `Z`, `initial` must contain two finite increasing numbers.",
        call. = FALSE
      )
    }
  } else if (length(initial) != length(z) + 1L) {
    stop(
      sprintf("With `Z`, `initial` must contain %d starting values.", length(z) + 1L),
      call. = FALSE
    )
  }
  as.numeric(initial)
}

.record_analysis_setup <- function(fit, call, z, stratum_setup) {
  fit$call <- call
  fit$details$z <- z
  fit$details$stratum <- if (stratum_setup$default) "constant" else stratum_setup$stratum
  fit
}

update_linear <- function(
    data, z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
    x = "X", x_star = "X_star", complete = "R", stratum = NULL,
    x_models = NULL, h = NULL, initial = NULL,
    n_quad = 30L, nboot = 0, seed = NULL,
    var_method = "standard") {
  x_models_missing <- missing(x_models)
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  .joint_validate_nboot(nboot)
  var_method <- .validate_var_method(var_method)
  internal_method <- if (var_method == "standard") "ordinary" else "fully_bootstrap"
  z <- .analysis_z_names(data, z)
  if (x_models_missing) x_models <- if (length(z)) "linear" else "kernel"
  x_models <- .analysis_x_models(z, x_models)
  initial <- if (length(x_models)) .analysis_initial(initial, z) else NULL
  stratum_setup <- .analysis_stratum(data, stratum)
  data <- stratum_setup$data
  if (length(z) == 0L) {
    fit <- .update_linear_optimal(
      data = data, y = y, x = x, x_star = x_star, complete = complete,
      strata = stratum_setup$stratum, beta_interval = initial, h = h,
      sdy = NULL, estimate_kernel = "kernel" %in% x_models,
      nboot = nboot, se_method = internal_method
    )
  } else {
    fit <- .update_linear_joint(
      data, z, y, x, x_star, complete, stratum_setup$stratum,
      x_models, initial, n_quad, nboot, internal_method
    )
  }
  fit <- .record_analysis_setup(fit, call, z, stratum_setup)
  fit$details$x_models <- x_models
  fit$details$var_method <- var_method
  fit$tune_parameters <- if (length(z) == 0L) {
    list(h = fit$details$h)
  } else {
    list()
  }
  fit
}

update_logistic <- function(
    data, z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
    x = "X", x_star = "X_star", complete = "R", stratum = NULL,
    x_models = NULL, h = NULL, initial = NULL,
    tune_para = 0.01, nboot_tune = 50L,
    n_quad = 30L, nboot = 0, seed = NULL,
    var_method = "standard") {
  x_models_missing <- missing(x_models)
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  .joint_validate_nboot(nboot)
  var_method <- .validate_var_method(var_method)
  internal_method <- if (var_method == "standard") "ordinary" else "fully_bootstrap"
  z <- .analysis_z_names(data, z)
  if (x_models_missing) x_models <- if (length(z)) "linear" else "kernel"
  x_models <- .analysis_x_models(z, x_models)
  initial <- if (length(x_models)) .analysis_initial(initial, z) else NULL
  stratum_setup <- .analysis_stratum(data, stratum)
  data <- stratum_setup$data
  tune_para <- .validate_tune_para(tune_para)

  if (length(z) == 0L) {
    selected_tune <- if (!("kernel" %in% x_models)) {
      NULL
    } else if (length(tune_para) == 1L) {
      tune_para[1L]
    } else {
      .select_logistic_tune_optimal(
        data, y, x, x_star, complete, stratum_setup$stratum,
        h, initial, tune_para, nboot_tune
      )
    }
    fit <- .update_logistic_optimal(
      data = data, y = y, x = x, x_star = x_star, complete = complete,
      strata = stratum_setup$stratum, beta_interval = initial, h = h,
      tune = selected_tune, estimate_kernel = "kernel" %in% x_models,
      nboot = nboot, se_method = internal_method
    )
  } else {
    selected_tune <- if (length(x_models) == 0L) {
      NULL
    } else if (length(tune_para) == 1L) {
      stats::setNames(rep(tune_para[1L], length(x_models)), x_models)
    } else {
      .select_logistic_tune_joint(
        data, z, y, x, x_star, complete, stratum_setup$stratum,
        x_models, tune_para, nboot_tune, n_quad, initial
      )
    }
    fit <- .update_logistic_joint(
      data, z, y, x, x_star, complete, stratum_setup$stratum,
      x_models, selected_tune, initial, n_quad, nboot, internal_method
    )
  }
  fit <- .record_analysis_setup(fit, call, z, stratum_setup)
  fit$details$x_models <- x_models
  fit$details$var_method <- var_method
  fit$details$tune_para <- selected_tune
  fit$tune_parameters <- if (length(z) == 0L) {
    list(h = fit$details$h, tune_para = selected_tune)
  } else {
    list(tune_para = selected_tune)
  }
  fit
}

update_cox <- function(
    data, z = grep("^Z[0-9]+$", names(data), value = TRUE), time = "Y",
    event = "D", x = "X", x_star = "X_star", complete = "R",
    stratum = NULL, x_models = NULL,
    h = NULL, initial = NULL, n_quad = 30L, nboot = 0, seed = NULL,
    var_method = "standard") {
  x_models_missing <- missing(x_models)
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  .joint_validate_nboot(nboot)
  var_method <- .validate_var_method(var_method)
  internal_method <- if (var_method == "standard") "ordinary" else "fully_bootstrap"
  z <- .analysis_z_names(data, z)
  if (x_models_missing) x_models <- if (length(z)) "linear" else "kernel"
  x_models <- .analysis_x_models(z, x_models)
  initial <- if (length(x_models)) .analysis_initial(initial, z) else NULL
  stratum_setup <- .analysis_stratum(data, stratum)
  data <- stratum_setup$data
  if (length(z) == 0L) {
    fit <- .update_cox_optimal(
      data = data, time = time, event = event, x = x, x_star = x_star,
      complete = complete, strata = stratum_setup$stratum,
      beta_interval = initial, h = h,
      estimate_kernel = "kernel" %in% x_models,
      nboot = nboot, se_method = internal_method
    )
  } else {
    fit <- .update_cox_joint(
      data, z, time, event, x, x_star, complete, stratum_setup$stratum,
      x_models, initial, n_quad, nboot, internal_method
    )
  }
  fit <- .record_analysis_setup(fit, call, z, stratum_setup)
  fit$details$x_models <- x_models
  fit$details$var_method <- var_method
  fit$tune_parameters <- if (length(z) == 0L) {
    list(h = fit$details$h)
  } else {
    list()
  }
  fit
}

.validate_tune_para <- function(tune_para) {
  if (!is.numeric(tune_para) || length(tune_para) < 1L || anyNA(tune_para) ||
      any(!is.finite(tune_para)) || any(tune_para < 0)) {
    stop("`tune_para` must contain non-negative finite numbers.", call. = FALSE)
  }
  sort(unique(as.numeric(tune_para)))
}

.validate_nboot_tune <- function(nboot_tune) {
  if (!is.numeric(nboot_tune) || length(nboot_tune) != 1L ||
      is.na(nboot_tune) || nboot_tune != as.integer(nboot_tune) ||
      nboot_tune < 2L) {
    stop("`nboot_tune` must be a whole number of at least 2.", call. = FALSE)
  }
  as.integer(nboot_tune)
}

.tuning_outlier_count <- function(draws) {
  draws <- as.matrix(draws)
  apply(draws, 2L, function(values) {
    center <- stats::median(values)
    scale <- stats::mad(values)
    if (!is.finite(scale) || scale == 0) {
      return(sum(abs(values - center) > sqrt(.Machine$double.eps)))
    }
    sum(abs((values - center) / scale) > 3.5)
  })
}

.collect_tuning_draws <- function(data, nboot_tune, fit_once, expected_length) {
  draws <- matrix(NA_real_, nrow = nboot_tune, ncol = expected_length)
  completed <- 0L
  attempts <- 0L
  max_attempts <- nboot_tune + max(10L, ceiling(nboot_tune / 4))
  while (completed < nboot_tune && attempts < max_attempts) {
    attempts <- attempts + 1L
    index <- sample.int(nrow(data), nrow(data), replace = TRUE)
    draw <- try(fit_once(data[index, , drop = FALSE]), silent = TRUE)
    if (!inherits(draw, "try-error") && length(draw) == expected_length &&
        all(is.finite(draw))) {
      completed <- completed + 1L
      draws[completed, ] <- as.numeric(draw)
    }
  }
  if (completed < nboot_tune) {
    stop(
      sprintf(
        "Only %d of %d tuning bootstrap fits succeeded.",
        completed, nboot_tune
      ),
      call. = FALSE
    )
  }
  draws
}

.select_logistic_tune_optimal <- function(
    data, y, x, x_star, complete, strata, h,
    initial, tune_para, nboot_tune) {
  nboot_tune <- .validate_nboot_tune(nboot_tune)
  for (candidate in tune_para) {
    fit_once <- function(boot_data) {
      fit <- .logistic_update_once(
        data = boot_data, y = y, x = x, x_star = x_star,
        complete = complete, strata = strata, beta_interval = initial,
        h = h, tune = candidate, estimate_kernel = TRUE
      )
      auxiliary <- fit$details$kernel_auxiliary_estimates
      auxiliary[["complete"]] - auxiliary[["full"]]
    }
    draws <- .collect_tuning_draws(data, nboot_tune, fit_once, 1L)
    if (all(.tuning_outlier_count(draws) == 0L)) return(candidate)
  }
  tune_para[length(tune_para)]
}

.joint_logistic_aux_delta <- function(
    data, z, y, x, x_star, complete, strata, x_models,
    tune, quadrature, initial) {
  .check_required_columns(data, c(y, x, x_star, complete, strata, z))
  complete_index <- data[[complete]] == 1
  data_compl <- data[complete_index, , drop = FALSE]
  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[complete_index]
  original_fit <- suppressWarnings(stats::glm(
    .joint_formula(y, c(x, z)), data = data_compl,
    family = stats::binomial(), weights = ipw_compl
  ))
  est_compl <- as.numeric(stats::coef(original_fit))
  if (any(!is.finite(est_compl))) {
    stop("The tuning bootstrap outcome model did not converge.", call. = FALSE)
  }
  solver_initial <- if (is.null(initial)) est_compl else initial
  out <- matrix(NA_real_, nrow = length(solver_initial), ncol = length(x_models),
                dimnames = list(c(x, z), x_models))
  for (model_type in x_models) {
    x_fit <- .joint_fit_x_model(
      data_compl, x, x_star, z, ipw_compl, model_type
    )
    score_full <- function(beta) logistic_score_modelx(
      beta, data, x_fit$fit, x_fit$se, quadrature$nodes,
      quadrature$weights, tune = tune, z = z
    )
    score_compl <- function(beta) logistic_score_modelx(
      beta, data_compl, x_fit$fit, x_fit$se, quadrature$nodes,
      quadrature$weights, ipw_compl, tune, z
    )
    est_full <- .joint_solve_score(solver_initial, score_full)
    est_compl <- .joint_solve_score(solver_initial, score_compl)
    out[, model_type] <- est_full - est_compl
  }
  out
}

.select_logistic_tune_joint <- function(
    data, z, y, x, x_star, complete, strata, x_models,
    tune_para, nboot_tune, n_quad, initial) {
  nboot_tune <- .validate_nboot_tune(nboot_tune)
  quadrature <- .joint_quadrature(n_quad)
  p <- length(z) + 1L
  counts <- array(
    NA_integer_, dim = c(p, length(x_models), length(tune_para)),
    dimnames = list(c(x, z), x_models, as.character(tune_para))
  )
  for (j in seq_along(tune_para)) {
    candidate <- tune_para[j]
    fit_once <- function(boot_data) {
      .joint_logistic_aux_delta(
        boot_data, z, y, x, x_star, complete, strata,
        x_models, candidate, quadrature, initial
      )
    }
    draws <- .collect_tuning_draws(
      data, nboot_tune, fit_once, p * length(x_models)
    )
    dim(draws) <- c(nboot_tune, p, length(x_models))
    for (k in seq_along(x_models)) {
      model_draws <- matrix(draws[, , k], nrow = nboot_tune, ncol = p)
      counts[, k, j] <- .tuning_outlier_count(model_draws)
    }
    if (all(counts[, , j] == 0L)) break
  }
  selected <- vapply(seq_along(x_models), function(k) {
    candidate_ok <- vapply(seq_along(tune_para), function(j) {
      values <- counts[, k, j]
      all(!is.na(values)) && all(values == 0L)
    }, logical(1))
    if (any(candidate_ok)) min(tune_para[candidate_ok]) else tune_para[length(tune_para)]
  }, numeric(1))
  stats::setNames(selected, x_models)
}
