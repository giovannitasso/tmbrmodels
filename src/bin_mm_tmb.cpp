// Keep includes minimal: TMB.hpp should handle dependencies.
#include <TMB.hpp>
#include <TMB.hpp> // Only include TMB

// Use the density namespace for MVNORM
using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
@@ -19,25 +21,46 @@ Type objective_function<Type>::operator() () {
  PARAMETER_VECTOR(betas);  // fixed effects
  PARAMETER_VECTOR(u);      // random effects

  // Parameters for unstructured covariance matrix
  // Parameters for manual covariance construction
  // Assumes 2 random effects per group (intercept and slope)
  PARAMETER_VECTOR(log_stdevs);   // log(stdev) vector [log(sigma_int), log(sigma_slope)]
  PARAMETER(transf_corr);      // Transformed correlation parameter (-inf, +inf)

  PARAMETER_VECTOR(log_stdevs); // log(stdev) vector [log(sigma_int), log(sigma_slope)]
  PARAMETER(transf_rho);      // Transformed correlation parameter (-inf, +inf), we'll use tanh

//==========================
// PRELIMINARY CALCULATIONS
//==========================

  // Setup for MVN density using TMB's built-in tools
  int n_reff_per_group = 2; // Hardcoded for intercept and slope for now
  density::UNSTRUCTURED_CORR_t<Type> nldens = density::UNSTRUCTURED_CORR_t<Type>(log_stdevs, transf_corr);
  // --- Manual Covariance Matrix Construction ---
  int n_reff_per_group = 2; // Hardcoded for intercept and slope

  // Linear predictor (using TMB/Eigen matrix operations)
  vector<Type> eta = X * betas + Z * u;
  // 1. Calculate standard deviations and correlation
  vector<Type> stdevs = exp(log_stdevs);
  Type rho = tanh(transf_rho); // rho is between -1 and 1
  
  // 2. Construct the Cholesky factor (L_corr) of the correlation matrix R = [[1, rho], [rho, 1]]
  matrix<Type> L_corr(n_reff_per_group, n_reff_per_group);
  L_corr(0,0) = 1.0;
  L_corr(0,1) = 0.0;
  L_corr(1,0) = rho;
  L_corr(1,1) = sqrt(1.0 - rho*rho);
  
  // 3. Construct the diagonal matrix S of standard deviations
  matrix<Type> S = matrix<Type>(n_reff_per_group, n_reff_per_group).setZero();
  S(0,0) = stdevs(0);
  S(1,1) = stdevs(1);

  // Apply the inverse-logit link function to get probabilities
  // Explicit cast to vector<Type> to avoid Eigen lazy-evaluation issues
  // 4. Construct the Cholesky factor (L_cov) of the covariance matrix Sigma = S * R * S^T = (S * L_corr) * (S * L_corr)^T
  matrix<Type> L_cov = S * L_corr;
  
  // 5. Create the MVN density object using the Cholesky factor L_cov
  MVNORM_t<Type> neg_log_dmvnorm(L_cov);
  // Alternative using full Sigma matrix (less efficient but maybe more stable?):
  // matrix<Type> Sigma = L_cov * L_cov.transpose();
  // MVNORM_t<Type> neg_log_dmvnorm(Sigma); 
  
  // --- Linear predictor and probability ---
  vector<Type> eta = X * betas + Z * u;
  // Explicit cast to avoid lazy evaluation issues
  vector<Type> p = vector<Type>(1.0 / (1.0 + exp(-eta))); 

//==========================
@@ -52,11 +75,11 @@ Type objective_function<Type>::operator() () {
    vector<Type> u_group_i(n_reff_per_group);
    u_group_i(0) = u(i);              // Intercept for group i
    u_group_i(1) = u(i + n_groups);   // Slope for group i
    nll += nldens(u_group_i); // Adds the *negative* log-density from the MVN prior
    // Use the MVN density object created above
    nll += neg_log_dmvnorm(u_group_i); // Adds the *negative* log-density
  }

  // Likelihood for the response: y_i ~ Binomial(1, p_i)
  // Use TMB's vectorized dbinom, subtracting log-prob because it's NLL
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
@@ -68,9 +91,6 @@ Type objective_function<Type>::operator() () {
  REPORT(u);

  // Report the standard deviations and correlation derived from parameters
  vector<Type> stdevs = exp(log_stdevs);
  Type rho = nldens.corr(0, 1); // Extract correlation between effect 1 and 2
  
  REPORT(stdevs);
  REPORT(rho);
