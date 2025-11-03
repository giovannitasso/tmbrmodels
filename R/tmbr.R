#' @title Fit a GLMM using a formula interface
#' @description This is the main user-facing function of the tmbrmodels package.
#' It allows fitting Generalized Linear Mixed Models (GLMMs) using a standard R formula.
#'
#' @param formula A two-sided linear formula object describing the fixed-effects and
#'   random-effects parts of the model. The structure is `response ~ fixed_effects + (random_effects | grouping_factor)`.
#'   Only one grouping factor is currently supported.
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
#'   formula = y ~ x + (x | group), # Intercept and slope random effect
#'   data = sim_data,
#'   family = "binomial"
#' )
#'
#' print(fit)
#' 
#' # --- Fit a model with only random intercepts ---
#' fit_intercept <- tmbr(
#'   formula = y ~ x + (1 | group), # Only intercept random effect
#'   data = sim_data,
#'   family = "binomial"
#' )
#' print(fit_intercept)
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
  
  # Extract group info and number of random effects per group
  if (length(lmod$reTrms$flist) > 1) {
      stop("Currently only one grouping factor is supported.")
  }
  n_groups <- nlevels(lmod$reTrms$flist[[1]])
  # n_reff_per_group is the number of columns in the model matrix for random effects per group
  # which corresponds to the number of terms in the random effects formula part (e.g., (1|g)=1, (x|g)=2)
  n_reff_per_group <- ncol(lmod$reTrms$Zt) / n_groups 
  
  # Get names for random effects from Zt column names if possible
  re_names <- colnames(lmod$reTrms$Zt)
  if (!is.null(re_names)) {
      # Extract unique names (remove group prefix if present)
      re_names <- unique(gsub(paste0(names(lmod$reTrms$flist)[1]), "", re_names, fixed = TRUE))
      # Clean potential leading ':' or empty strings if intercept was implicit
      re_names <- sub("^:", "", re_names)
      if (any(re_names == "")) re_names[re_names == ""] <- "(Intercept)"
  } else {
      re_names <- NULL # Will use defaults later if names are missing
  }


  # --- 3. Call the model fitting engine ---
  fit <- fit_binomial_glmm(y = y, 
                           X = X, 
                           Z = Z, 
                           n_groups = n_groups, 
                           n_reff_per_group = n_reff_per_group,
                           re_names = re_names) # Pass names down

  # --- 4. Return the result ---
  return(fit)
}
