.update_cox_joint <- function(data, z, time, event, x, x_star, complete,
                              strata, x_models, initial, n_quad, nboot,
                              se_method = "ordinary") {
  call <- match.call()
  .joint_validate_nboot(nboot)
  z <- .joint_z_names(data, z)
  x_models <- .joint_x_models(x_models)
  quadrature <- .joint_quadrature(n_quad)

  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .cox_joint_once(
    data, z, time, event, x, x_star, complete, strata, x_models,
    quadrature, initial, compute_update = !fully_bootstrap
  )
  fit_once <- function(boot_data) {
    .cox_joint_once(
      boot_data, z, time, event, x, x_star, complete, strata, x_models,
      quadrature, initial, compute_update = !fully_bootstrap
    )[[if (fully_bootstrap) "components" else "estimates"]]
  }
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot, fit_once, fit$components
    )
    fit$estimates <- bootstrap_fit$estimates
    se <- bootstrap_fit$se
    fit$details$projection <- stats::setNames(
      rep("full_bootstrap", ncol(fit$estimates) - 1L),
      colnames(fit$estimates)[-1L]
    )
    fit$details$full_bootstrap_covariance <- bootstrap_fit$covariance
    fit$details$bootstrap_attempts <- bootstrap_fit$attempts
    fit$details$bootstrap_failures <- bootstrap_fit$failures
  } else {
    se <- .joint_bootstrap_se(data, nboot, fit_once, fit$estimates)
  }
  .new_joint_fit(call, "cox", fit$estimates, se, fit$details)
}

.cox_joint_once <- function(data, z, time, event, x, x_star, complete,
                            strata, x_models, quadrature, initial,
                            compute_update = TRUE) {
  .check_required_columns(data, c(time, event, x, x_star, complete, strata, z))
  complete_index <- data[[complete]] == 1
  data_compl <- data[complete_index, , drop = FALSE]
  p <- length(z) + 1L
  if (nrow(data_compl) <= p || sum(data_compl[[event]] == 1) == 0) {
    stop("The Phase II sample is too small for the joint Cox model.", call. = FALSE)
  }
  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[complete_index]
  centered <- data[[complete]] * ipw_all - 1
  terms <- c(x, z)
  response <- paste0("survival::Surv(", time, ", ", event, ")")
  original_formula <- .joint_formula(response, terms)
  working_formula <- .joint_formula(response, c(x_star, z))

  original_fit <- survival::coxph(original_formula, data = data_compl,
                                  weights = ipw_compl)
  working_fit <- survival::coxph(working_formula, data = data)
  working_fit_compl <- survival::coxph(working_formula, data = data_compl,
                                       weights = ipw_compl)
  est_compl <- as.numeric(stats::coef(original_fit))
  est_work <- as.numeric(stats::coef(working_fit))
  est_work_compl <- as.numeric(stats::coef(working_fit_compl))
  if (any(!is.finite(c(est_compl, est_work, est_work_compl)))) {
    stop("A joint Cox model is rank deficient or did not converge.", call. = FALSE)
  }
  solver_initial <- if (is.null(initial)) est_compl else initial

  c_compl <- .joint_design(data_compl, x, z)
  differences <- list(default = est_work_compl - est_work)
  methods <- list(complete_case = est_compl)
  projection_methods <- character()

  if (compute_update) {
    info_compl <- .cox_joint_info_org(
      est_compl, data, data_compl, ipw_compl, time, event, x, z
    )
    info_work <- .cox_joint_info_work(est_work, data, time, event, x_star, z)
    score_compl <- .cox_joint_score_org_vec(
      est_compl, data_compl, data, data_compl, ipw_compl,
      rep(1, nrow(data_compl)), time, event, x, z
    )
    raw_compl <- matrix(0, nrow(data), p)
    raw_compl[complete_index, ] <- score_compl * ipw_compl
    score_work <- .cox_joint_score_work_vec(
      est_work, data, data, rep(1, nrow(data)), time, event, x_star, z
    )
    influence_compl <- .joint_influence(raw_compl, info_compl)
    influence_work <- .joint_influence(score_work * centered, info_work)
    outcome_projection <- .joint_projection(influence_work, influence_compl)
    methods$default <- as.numeric(
      est_compl - outcome_projection$coef %*% differences$default
    )
    projection_methods["default"] <- outcome_projection$method
  }

  add_conditional_update <- function(label, vector_full, vector_compl) {
    score_full <- function(beta) colSums(vector_full(beta))
    score_phase2 <- function(beta) colSums(vector_compl(beta))
    est_full <- .joint_solve_score(solver_initial, score_full)
    est_phase2 <- .joint_solve_score(solver_initial, score_phase2)
    difference <- est_phase2 - est_full
    joint_label <- paste0("joint_", label)
    differences[[label]] <<- difference
    differences[[joint_label]] <<- c(differences$default, difference)
    if (!compute_update) return(invisible(NULL))

    information <- -numDeriv::jacobian(score_full, est_full)
    influence <- .joint_influence(vector_full(est_full) * centered, information)
    single_projection <- .joint_projection(influence, influence_compl)
    joint_projection <- .joint_projection(
      cbind(influence_work, influence), influence_compl
    )
    methods[[label]] <<- as.numeric(
      est_compl - single_projection$coef %*% difference
    )
    methods[[joint_label]] <<- as.numeric(
      est_compl - joint_projection$coef %*%
        differences[[joint_label]]
    )
    projection_methods[label] <<- single_projection$method
    projection_methods[joint_label] <<- joint_projection$method
  }

  for (model_type in x_models) {
    x_fit <- .joint_fit_x_model(data_compl, x, x_star, z, ipw_compl, model_type)
    pred_all <- as.numeric(stats::predict(x_fit$fit, data))
    pred_compl <- as.numeric(stats::predict(x_fit$fit, data_compl))
    vector_full <- function(beta) .cox_joint_conditional_score_vec(
      beta, data, data, data_compl, c_compl, ipw_compl,
      rep(1, nrow(data)), quadrature$nodes, quadrature$weights,
      time, event, z, pred_x = pred_all,
      se_x = x_fit$se, est_compl = est_compl
    )
    vector_compl <- function(beta) .cox_joint_conditional_score_vec(
      beta, data_compl, data, data_compl, c_compl, ipw_compl,
      ipw_compl, quadrature$nodes, quadrature$weights,
      time, event, z, pred_x = pred_compl,
      se_x = x_fit$se, est_compl = est_compl
    )
    add_conditional_update(model_type, vector_full, vector_compl)
  }

  components <- .bootstrap_components(
    stats::setNames(est_compl, terms), differences
  )
  estimates <- do.call(cbind, methods)
  rownames(estimates) <- terms
  list(
    estimates = estimates,
    components = components,
    details = list(
      z = z, ipw = ipw_all, projection = projection_methods,
      x_models = x_models, initial = initial,
      quadrature_points = length(quadrature$nodes)
    )
  )
}
