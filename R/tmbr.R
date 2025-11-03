#' @title Fit a GLMM using a formula interface
#' @description This is the main user-facing function of the tmbrmodels package.
#' It allows fitting Generalized Linear Mixed Models (GLMMs) using a standard R formula.
#'
#' @param formula A two-sided linear formula object describing the fixed-effects and
#'   random-effects parts of the model. The structure is `response ~ fixed_effects + (random_effects | grouping_factor)`.
#' @param data A data frame containing the variables named in the formula.
#' @param family Currently, only `family = "binomial"` is supported.
#'
#' @return An object of class `tmbr_fit`.
#' @import Matrix
#' @importFrom stats model.response
#' @importFrom methods as
#' @export
#' @examples
#' \dontrun{
#' # Simulate data for the formula interface
#' set.seed(1234)
#' n_groups <- 10
#' n_per_group <- 30
#'
#' # Create a data frame, which is required for formula interfaces
#' sim_data <- data.frame(
#'   group = factor(rep(1:n_groups, each = n_per_group)),
#'   x = rnorm(n_groups * n_per_group)
#' )
#'
#' # Simulate random effects for intercepts and slopes
#' u_intercepts <- rnorm(n_groups, 0, 1.0)
#' u_slopes <- rnorm(n_groups, 0, 0.5)
#'
#' # Calculate the linear predictor to generate the response variable
#' # Fixed effects: intercept=-1.0, slope=1.5
#' eta <- (-1.0 + u_intercepts[sim_data$group]) +
#'        (1.5 + u_slopes[sim_data$group]) * sim_data$x
#'
#' p <- 1 / (1 + exp(-eta))
#' sim_data$y <- rbinom(nrow(sim_data), size = 1, prob = p)
#'
#' # --- Fit the model using the new, simple formula interface! ---
#' library(tmbrmodels)
#' fit <- tmbr(
#'   formula = y ~ x + (x | group),
#'   data = sim_data,
#'   family = "binomial"
#' )
#'
#' print(fit)
#' }
tmbr <- function(formula, data, family = "binomial") {

  # --- 1. Initial Validations ---
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("The 'formula' argument must be a valid R formula.")
  }
  if (missing(data)) {
    stop("The 'data' argument (a data frame) must be provided.")
  }
  if (family != "binomial") {
    stop("Currently, only the 'binomial' family is supported.")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The 'lme4' package is required. Please install it with: install.packages('lme4')")
  }

  # --- 2. Formula Parsing with lme4 ---
  lmod <- lme4::lFormula(formula = formula, data = data)

  y <- stats::model.response(lmod$fr)
  X <- lmod$X
  Z <- as(t(lmod$reTrms$Zt), "matrix")
  
  # Extract the number of groups
  n_groups <- nlevels(lmod$reTrms$flist[[1]])

  # --- 3. Call the model fitting engine ---
  fit <- fit_binomial_glmm(y = y, X = X, Z = Z, n_groups = n_groups)

  # --- 4. Return the result ---
  return(fit)
}
