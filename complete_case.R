#' Complete-Case analysis
#' 
#' @param y
#' @param X, model matrix including a column for intercept
#' @return list of point estimates and standard errors
complete_case_lm <- function(y, X) {
  p <- ncol(X) - 1
  n <- nrow(na.omit(X))
  # intercept include in X
  model <- lm(y ~ 0 + X)
  # names
  var_names <- c(paste("mu", seq_len(p), sep = "_"),
                 paste("sigma", seq_len(p), "sq", sep = "_"),
                 "sigma_sq",
                 paste("b", seq(0,p), sep = "_"))
  # point estimation
  mu <- unname(apply(X, 2, mean, na.rm = TRUE)[-1])
  Sigma_x <- unname(var(X, na.rm = TRUE)[-1,-1])
  sigma_sq <- unname(mean(residuals(model)^2))
  b <- unname(coef(model))
  point_est <- list(mu = mu, Sigma_x = Sigma_x, sigma_sq = sigma_sq, b = b)
  # standard error
  mu.se <- unname(apply(X, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))[-1])
  Sigma_x.diag.se <- unname(apply(X, 2, function(x) var(x, na.rm = TRUE) * sqrt(2/(length(na.omit(x)) - 1)))[-1])
  sigma_sq.se <- sigma_sq * sqrt(2 * (n - p) / n ^ 2)
  b.se <- unname(sqrt(diag(vcov(model))))
  se <- c(mu.se, Sigma_x.diag.se, sigma_sq.se, b.se)
  names(se) <- var_names
  
  return(list(point_est = point_est, se = se))
}
