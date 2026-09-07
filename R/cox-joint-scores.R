.cox_joint_named_matrix <- function(data, time, event, x = NULL,
                                    x_star = NULL, z) {
  n <- nrow(data)
  x_value <- if (is.null(x)) numeric(n) else data[[x]]
  x_star_value <- if (is.null(x_star)) numeric(n) else data[[x_star]]
  out <- cbind(
    complete = numeric(n), time = data[[time]], event = data[[event]],
    x = x_value, x_star = x_star_value,
    as.matrix(data[, z, drop = FALSE])
  )
  storage.mode(out) <- "double"
  out
}

.cox_joint_info_org <- function(beta, data_all, data_compl, ipw_compl,
                                time, event, x, z) {
  cpp_cox_joint_info_org(
    beta,
    .cox_joint_named_matrix(data_all, time, event, z = z),
    .cox_joint_named_matrix(data_compl, time, event, x = x, z = z),
    ipw_compl
  )
}

.cox_joint_info_work <- function(beta, data_all, time, event, x_star, z) {
  cpp_cox_joint_info_work_coxph(
    beta,
    .cox_joint_named_matrix(data_all, time, event, x_star = x_star, z = z)
  )
}

.cox_joint_score_org_vec <- function(beta, data, data_all, data_compl,
                                     ipw_compl, ipw, time, event, x, z) {
  cpp_cox_joint_score_org_vec(
    beta,
    .cox_joint_named_matrix(data, time, event, x = x, z = z),
    .cox_joint_named_matrix(data_all, time, event, z = z),
    .cox_joint_named_matrix(data_compl, time, event, x = x, z = z),
    ipw_compl, ipw
  )
}

.cox_joint_score_work_vec <- function(beta, data, data_all, ipw,
                                      time, event, x_star, z) {
  cpp_cox_joint_score_work_coxph_vec(
    beta,
    .cox_joint_named_matrix(data, time, event, x_star = x_star, z = z),
    .cox_joint_named_matrix(data_all, time, event, x_star = x_star, z = z),
    ipw
  )
}

.cox_joint_conditional_score_vec <- function(beta, data, data_all, data_compl,
                                             complete_covariates, ipw_compl, ipw,
                                             grid, weight, time, event, z,
                                             pred_x, se_x, est_compl) {
  n_compl <- nrow(data_compl)
  complete_matrix <- cbind(
    complete = numeric(n_compl), time = data_compl[[time]],
    event = data_compl[[event]], x = complete_covariates[, 1L],
    x_star = numeric(n_compl), complete_covariates[, -1L, drop = FALSE]
  )
  storage.mode(complete_matrix) <- "double"
  cpp_cox_joint_score_work_modelX_vec(
    beta,
    .cox_joint_named_matrix(data, time, event, z = z),
    .cox_joint_named_matrix(data_all, time, event, z = z),
    complete_matrix, pred_x, se_x, grid, weight, pi,
    est_compl, ipw_compl, ipw
  )
}
