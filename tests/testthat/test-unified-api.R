test_that("the package exposes one generator and updater per family", {
  expect_identical(packageName(asNamespace("twoStep2Phase")), "twoStep2Phase")
  expect_setequal(
    getNamespaceExports("twoStep2Phase"),
    c(
      "generate_linear_data", "generate_logistic_data", "generate_cox_data",
      "update_linear", "update_logistic", "update_cox"
    )
  )
})

test_that("the analysis seed is set once at updater entry", {
  for (fun in list(update_linear, update_logistic, update_cox)) {
    body_text <- paste(deparse(body(fun)), collapse = "\n")
    expect_equal(lengths(regmatches(
      body_text, gregexpr("set.seed(seed)", body_text, fixed = TRUE)
    )), 1L)
  }
  expect_false(grepl(
    "seed = NULL", paste(deparse(body(update_logistic)), collapse = "\n"),
    fixed = TRUE
  ))

  internal_functions <- c(
    ".update_linear_optimal", ".update_logistic_optimal", ".update_cox_optimal",
    ".update_linear_joint", ".update_logistic_joint", ".update_cox_joint",
    ".joint_bootstrap_se", ".fully_bootstrap_fit"
  )
  for (name in internal_functions) {
    fun <- getFromNamespace(name, "twoStep2Phase")
    expect_false("seed" %in% names(formals(fun)))
    expect_false(grepl("set.seed", paste(deparse(body(fun)), collapse = "\n"),
                       fixed = TRUE))
  }
})

test_that("linear quantiles define low-to-high sampling strata", {
  set.seed(10)
  data <- generate_linear_data(
    100, 30, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1,
    ns = c(6, 15, 9), quantile = c(0.3, 0.7)
  )

  expect_equal(sort(unique(data$stratum)), 0:2)
  selected <- table(factor(data$stratum[data$R == 1], levels = 0:2))
  expect_equal(as.numeric(selected), c(6, 15, 9))
  expect_equal(sum(data$R), 30)
  expect_true(all(is.na(data$X[data$R == 0])))
})

test_that("linear quantile NULL uses one sampling stratum", {
  set.seed(11)
  no_z <- generate_linear_data(
    100, 30, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1
  )
  with_z <- generate_linear_data(
    100, 30, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 1
  )

  expect_equal(unique(no_z$stratum), 0)
  expect_equal(sum(no_z$R), 30)
  expect_false("Z1" %in% names(no_z))
  expect_true("Z1" %in% names(with_z))
  expect_false(anyNA(with_z$X_star))
})

test_that("generator inputs are validated", {
  expect_error(
    generate_linear_data(
      100, 30, beta = c(0.5, 0.25), mu = 0, Sigma = matrix(1),
      sde = 1
    ),
    "length 2"
  )
  expect_error(
    generate_linear_data(
      100, 30, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1,
      ns = c(10, 20), quantile = c(0.3, 0.7)
    ),
    "one sample size per stratum"
  )
})

test_that("analysis defaults depend on whether Z is present", {
  set.seed(12)
  no_z <- generate_linear_data(
    100, 30, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  no_z$stratum <- NULL
  fit_no_z <- update_linear(no_z)

  with_z <- generate_linear_data(
    100, 30, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )
  with_z$stratum <- NULL
  fit_with_z <- update_linear(with_z, n_quad = 5, seed = 1)

  expect_identical(names(fit_no_z$estimate), c("complete_case", "default", "kernel"))
  expect_identical(colnames(fit_with_z$estimate),
                   c("complete_case", "default", "linear", "joint_linear"))
  expect_identical(fit_no_z$details$stratum, "constant")
  expect_identical(fit_with_z$details$z, "Z1")
  expect_identical(fit_with_z$details$quadrature_points, 5L)
  expect_false("n_quad" %in% names(fit_with_z$tune_parameters))
  expect_identical(names(fit_no_z$tune_parameters), "h")
  expect_identical(fit_with_z$tune_parameters, list())
})

test_that("incompatible x_models produce errors and NULL omits them", {
  set.seed(13)
  no_z <- generate_linear_data(
    100, 30, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  with_z <- generate_linear_data(
    100, 30, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )

  expect_error(update_linear(no_z, x_models = "linear"), "Without `Z`")
  expect_error(update_linear(with_z, x_models = "kernel"), "With `Z`")
  expect_identical(
    names(update_linear(no_z, x_models = NULL)$estimate),
    c("complete_case", "default")
  )
  expect_identical(
    colnames(update_linear(with_z, x_models = NULL, seed = 1)$estimate),
    c("complete_case", "default")
  )
})

test_that("explicit NULL forces analysis without Z", {
  set.seed(131)
  data <- generate_linear_data(
    120, 40, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )

  fit <- update_linear(data, z = NULL, x_models = NULL)
  expect_identical(fit$details$z, character())
  expect_identical(names(fit$estimate), c("complete_case", "default"))
  expect_identical(
    twoStep2Phase:::.analysis_z_names(data, character(0)),
    character()
  )

  fit_null <- update_linear(data, z = NULL)
  fit_empty <- update_linear(data, z = character(0))
  expect_identical(names(fit_null$estimate),
                   c("complete_case", "default", "kernel"))
  expect_identical(names(fit_empty$estimate),
                   c("complete_case", "default", "kernel"))
})

test_that("removed update arguments are absent", {
  for (fun in list(update_linear, update_logistic, update_cox)) {
    expect_false("beta_interval" %in% names(formals(fun)))
    expect_false("estimate_kernel" %in% names(formals(fun)))
  }
  expect_false("sdy" %in% names(formals(update_linear)))
  for (fun in list(update_linear, update_logistic, update_cox)) {
    expect_true("initial" %in% names(formals(fun)))
    expect_false("intv" %in% names(formals(fun)))
    expect_true("stratum" %in% names(formals(fun)))
    expect_false("strata" %in% names(formals(fun)))
    expect_true("var_method" %in% names(formals(fun)))
    expect_false("se_method" %in% names(formals(fun)))
    expect_identical(formals(fun)$var_method, "standard")
  }
  expect_identical(formals(update_linear)$complete, "R")
})

test_that("initial accepts bounds without Z and starting values with Z", {
  set.seed(132)
  no_z <- generate_linear_data(
    120, 50, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  expect_true(all(is.finite(update_linear(
    no_z, z = NULL, initial = c(0.1, 0.9), nboot = 0
  )$estimate)))

  joint_data <- generate_linear_data(
    140, 70, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )
  expect_true(all(is.finite(update_linear(
    joint_data, initial = c(0.5, 0.25), n_quad = 5, nboot = 0
  )$estimate)))
  expect_error(update_linear(joint_data, initial = 0.5), "2 starting values")

  set.seed(133)
  logistic_data <- generate_logistic_data(
    180, 110, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )
  logistic_fit <- suppressWarnings(update_logistic(
    logistic_data, initial = c(0.5, 0.25), n_quad = 5, nboot = 0,
    tune_para = 0.01
  ))
  expect_true(all(is.finite(logistic_fit$estimate)))
  expect_equal(logistic_fit$details$initial, c(0.5, 0.25))

  set.seed(134)
  cox_data <- generate_cox_data(
    180, 110, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )
  cox_fit <- update_cox(
    cox_data, initial = c(0.5, 0.25), n_quad = 5, nboot = 0
  )
  expect_true(all(is.finite(cox_fit$estimate)))
  expect_equal(cox_fit$details$initial, c(0.5, 0.25))
})

test_that("a singleton logistic tuning value skips tuning bootstrap", {
  set.seed(14)
  data <- generate_logistic_data(
    120, 30, beta = 0.5, mu = 0, Sigma = matrix(1),
    sde = 0.5, ns = c(15, 15)
  )
  data$stratum <- NULL

  fit <- update_logistic(data, tune_para = 0.01, nboot_tune = 1)
  expect_equal(fit$tune_parameters$tune_para, 0.01)
  expect_true(is.finite(fit$tune_parameters$h))
})

test_that("generators use sde and logistic and Cox default to MCAR", {
  expect_true("sde" %in% names(formals(generate_linear_data)))
  expect_false("error_var" %in% names(formals(generate_linear_data)))
  expect_identical(formals(generate_logistic_data)$ns, quote(n))
  expect_identical(formals(generate_cox_data)$ns, quote(n))
  expect_equal(formals(generate_cox_data)$tau, 2.5)

  set.seed(15)
  logistic <- generate_logistic_data(
    500, 150, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1
  )
  cox <- generate_cox_data(
    500, 150, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1
  )
  expect_identical(unique(logistic$stratum), 0L)
  expect_identical(unique(cox$stratum), 0L)
  expect_equal(sum(logistic$R), 150)
  expect_equal(sum(cox$R), 150)
})

test_that("sde is the measurement-error standard deviation", {
  set.seed(16)
  data <- generate_linear_data(
    5000, 5000, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 2
  )
  expect_equal(stats::sd(data$X_star - data$X), 2, tolerance = 0.08)
})

test_that("GAM working-model scale uses raw response residuals", {
  set.seed(17)
  data <- generate_linear_data(
    200, 100, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 1
  )
  complete <- data$R == 1
  weights_all <- twoStep2Phase:::.ipw_by_stratum(
    data, complete = "R", strata = "stratum"
  )
  weights <- weights_all[complete]
  fit <- twoStep2Phase:::.joint_fit_x_model(
    data[complete, , drop = FALSE], "X", "X_star", "Z1", weights, "gam"
  )
  expected <- sqrt(sum(weights * fit$fit$residuals^2) / sum(weights))

  expect_equal(fit$se, expected)
  expect_false(isTRUE(all.equal(
    fit$se,
    sqrt(sum(weights * stats::residuals(fit$fit)^2) / sum(weights))
  )))
})
