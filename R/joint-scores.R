.d_log_normal_pdf <- function(t, m, s) -(t - m) / (s^2)
.dd_log_normal_pdf <- function(s) -1 / (s^2)
.dm_log_normal_pdf <- function(t, m, s) (t - m) / (s^2)
.ddm_log_normal_pdf <- function(s) -1 / (s^2)

linear_score_modelx <- function(beta, data, resX, seX, sdy, grid, weight,
                                ipw = NULL, z = NULL) {
  .linear_score_modelx_impl(
    beta, data, resX, seX, sdy, grid, weight, ipw,
    vector = FALSE, z = z
  )
}

logistic_score_modelx <- function(beta, data, resX, seX, grid, weight, ipw = NULL,
                                  tune = NULL, z = NULL) {
  res <- .logistic_score_modelx_impl(beta, data, resX, seX, grid, weight, ipw, z)
  if (is.null(tune)) {
    colSums(res)
  } else {
    colMeans(res) - 2 * beta * (nrow(data)^(-1 / 3)) * tune
  }
}

.z_matrix <- function(data, z = NULL) {
  if (is.null(z)) z <- grep("^Z[0-9]+$", names(data), value = TRUE)
  as.matrix(data[, z, drop = FALSE])
}

.linear_score_modelx_impl <- function(beta, data, resX, seX, sdy, grid, weight,
                                      ipw, vector, z = NULL) {
  res <- matrix(0, nrow = nrow(data), ncol = length(beta))
  for (i in seq_len(nrow(data))) {
    Zi <- as.numeric(.z_matrix(data[i, , drop = FALSE], z))
    predXi <- stats::predict(resX, data[i, , drop = FALSE])
    b <- 0
    for (j in seq_len(100)) {
      old_b <- b
      score <- .dm_log_normal_pdf(data$Y[i], sum(c(b, Zi) * beta), sdy) * beta[1] +
        .d_log_normal_pdf(b, predXi, seX)
      hes <- .ddm_log_normal_pdf(sdy) * (beta[1]^2) + .dd_log_normal_pdf(seX)
      b <- b - score / hes
      if (abs(b - old_b) < 1e-3) break
    }
    gridb <- grid / sqrt(-hes) + b
    cova <- cbind(gridb, matrix(rep(Zi, length(gridb)), byrow = TRUE, nrow = length(gridb)))
    score_org <- cova * (data$Y[i] - as.vector(cova %*% beta))
    log_fY_XZ <- stats::dnorm(data$Y[i], as.vector(cova %*% beta), sdy, log = TRUE)
    log_fX_XstarZ <- stats::dnorm(gridb, predXi, seX, log = TRUE)
    weight1 <- exp(log_fY_XZ - max(log_fY_XZ) + log_fX_XstarZ - max(log_fX_XstarZ) +
                     grid^2 + log(weight))
    res[i, ] <- colSums(score_org * weight1) / sum(weight1)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  if (vector) res else colSums(res)
}

.logistic_score_modelx_impl <- function(beta, data, resX, seX, grid, weight,
                                        ipw, z = NULL) {
  res <- matrix(0, nrow = nrow(data), ncol = length(beta))
  for (i in seq_len(nrow(data))) {
    Zi <- as.numeric(.z_matrix(data[i, , drop = FALSE], z))
    predXi <- stats::predict(resX, data[i, , drop = FALSE])
    cova <- cbind(grid, matrix(rep(Zi, length(grid)), byrow = TRUE, nrow = length(grid)))
    eta <- as.vector(cova %*% beta)
    score_org <- cova * (data$Y[i] - stats::plogis(eta))
    log_fY_XZ <- eta * data$Y[i] - .log1pexp(eta)
    log_fX_XstarZ <- stats::dnorm(grid, predXi, seX, log = TRUE)
    weight1 <- exp(log_fY_XZ - max(log_fY_XZ) + log_fX_XstarZ - max(log_fX_XstarZ) +
                     grid^2 + log(weight))
    res[i, ] <- colSums(score_org * weight1) / sum(weight1)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  res
}
