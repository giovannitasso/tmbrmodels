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
  DATA_INTEGER(n_groups);         // Number of groups
  DATA_INTEGER(n_reff_per_group); // Number of random effects per group

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);        // Fixed effects
  PARAMETER_VECTOR(u);            // Random effects vector (all groups combined)
  PARAMETER_VECTOR(log_stdevs);   // Vector of log standard deviations for random effects
  PARAMETER_CHOLESKY_CORR(L_corr); // Cholesky factor of the correlation matrix (parameters handled by TMB)


//==========================
// PRELIMINARY CALCULATIONS
//==========================
  
  // --- Construct Covariance Matrix Components ---
  
  // 1. Cholesky factor of Correlation matrix R (comes directly from PARAMETER_CHOLESKY_CORR)
  // L_corr is lower triangular matrix of size n_reff_per_group x n_reff_per_group
  
  // 2. Diagonal matrix S of standard deviations
  vector<Type> stdevs = exp(log_stdevs);
  matrix<Type> S = matrix<Type>(n_reff_per_group, n_reff_per_group).setZero();
  for(int i = 0; i < n_reff_per_group; ++i) {
      S(i,i) = stdevs(i);
  }
  
  // 3. Cholesky factor (L_cov) of the covariance matrix Sigma = S * R * S^T = (S * L_corr) * (S * L_corr)^T
  matrix<Type> L_cov = S * L_corr;
  
  // 4. Create the MVN density object using the Cholesky factor L_cov
  MVNORM_t<Type> neg_log_dmvnorm(L_cov);
  
  // --- Linear predictor and probability ---
  vector<Type> eta = X * betas + Z * u;
  // Explicit cast to avoid lazy evaluation issues
  vector<Type> p = vector<Type>(1.0 / (1.0 + exp(-eta))); 

//==========================
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0; // Initialize negative log-likelihood

  // Prior for random effects: u ~ MVN(0, Sigma)
  // Loop through each random effect type, then each group
  vector<Type> u_group_i(n_reff_per_group); // Reusable vector for each group's effects
  for(int i = 0; i < n_groups; ++i){
      for(int j = 0; j < n_reff_per_group; ++j){
          // Assumes u is ordered [eff1_g1..N, eff2_g1..N, ..., effk_g1..N]
          u_group_i(j) = u(j * n_groups + i); 
      }
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
  
  // Report standard deviations and the full correlation matrix R
  matrix<Type> R = L_corr * L_corr.transpose();
  
  REPORT(stdevs);
  REPORT(R); // Report the correlation matrix
  
  // Calculate standard errors for derived quantities using the delta method
  ADREPORT(stdevs);
  ADREPORT(R);

  return nll;
}
