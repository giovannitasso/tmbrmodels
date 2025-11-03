#include <TMB.hpp> // Only include TMB

// Use the density namespace for MVNORM
using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
  
//==========================
// DATA SECTION
//==========================
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_INTEGER(n_groups); // Number of groups

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);  // fixed effects
  PARAMETER_VECTOR(u);      // random effects
  
  // Parameters for manual covariance construction
  // Assumes 2 random effects per group (intercept and slope)
  PARAMETER_VECTOR(log_stdevs); // log(stdev) vector [log(sigma_int), log(sigma_slope)]
  PARAMETER(transf_rho);      // Transformed correlation parameter (-inf, +inf), we'll use tanh

//==========================
// PRELIMINARY CALCULATIONS
//==========================
  
  // --- Manual Covariance Matrix Construction ---
  int n_reff_per_group = 2; // Hardcoded for intercept and slope
  
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
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0; // Initialize negative log-likelihood

  // Prior for random effects: u ~ MVN(0, Sigma)
  // Loop through each group and apply the MVN density
  // Assumes u is ordered as [int_g1, ..., int_gN, slope_g1, ..., slope_gN]
  for(int i = 0; i < n_groups; ++i){
    vector<Type> u_group_i(n_reff_per_group);
    u_group_i(0) = u(i);              // Intercept for group i
    u_group_i(1) = u(i + n_groups);   // Slope for group i
    // Use the MVN density object created above
    nll += neg_log_dmvnorm(u_group_i); // Adds the *negative* log-density
  }

  // Likelihood for the response: y_i ~ Binomial(1, p_i)
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  REPORT(betas);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  
  // Report the standard deviations and correlation derived from parameters
  REPORT(stdevs);
  REPORT(rho);
  
  // Calculate standard errors for derived quantities using the delta method
  ADREPORT(stdevs);
  ADREPORT(rho);

  return nll;
}
