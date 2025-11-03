#include <TMB.hpp>

// Define M_PI (pi) if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template<class Type>
Type objective_function<Type>::operator() () {
  
//==========================
// DATA SECTION
//==========================
  DATA_VECTOR(y);             // Response vector (0/1)
  DATA_MATRIX(X);             // Fixed effects design matrix
  DATA_MATRIX(Z);             // Random effects design matrix
  DATA_INTEGER(n_groups);         // Number of groups (e.g., subjects, sites)
  DATA_INTEGER(n_reff_per_group); // Number of random effects per group (e.g., 2 for intercept + slope)

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);        // Fixed effects coefficients
  PARAMETER_VECTOR(u);            // Random effects vector (all groups combined)
  PARAMETER_VECTOR(log_stdevs);   // Vector of log-standard deviations for random effects
  PARAMETER_CHOLESKY_CORR(L_corr); // Cholesky factor of the correlation matrix

//==========================
// PRELIMINARY CALCULATIONS
//==========================
  
  // --- Manual Covariance Construction ---
  
  // 1. Standard deviations
  vector<Type> stdevs = exp(log_stdevs);
  
  // 2. Diagonal matrix S of standard deviations
  matrix<Type> S = matrix<Type>(n_reff_per_group, n_reff_per_group).setZero();
  for(int i = 0; i < n_reff_per_group; ++i) {
      S(i,i) = stdevs(i);
  }
  
  // 3. Cholesky factor (L_cov) of the covariance matrix
  // Sigma = S * R * S^T = (S * L_corr) * (S * L_corr)^T
  // L_cov = S * L_corr
  matrix<Type> L_cov = S * L_corr;

  // 4. log-determinant of L_cov
  // log(det(L_cov)) = sum(log(diag(L_cov)))
  Type logdet_L_cov = 0.0;
  for (int i=0; i < n_reff_per_group; i++) {
      logdet_L_cov += log(L_cov(i,i));
  }
  
  // 5. Constant for the NLL
  Type log_2pi = log(2.0 * M_PI);

  // --- Linear Predictor and Probability ---
  vector<Type> eta = X * betas + Z * u;
  // Explicit cast to vector<Type> to avoid Eigen lazy-evaluation issues
  vector<Type> p = vector<Type>(1.0 / (1.0 + exp(-eta))); 

//==========================
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0; // NLL = Negative Log-Likelihood

  // --- Manual Prior for Random Effects ---
  // NLL = -log(L(u | Sigma))
  // NLL = k/2 * log(2pi) + log(det(L_cov)) + 0.5 * ||v||^2
  // where v = L_cov^{-1} * u
  
  vector<Type> u_group_i(n_reff_per_group); // Reusable vector for this group's effects
  vector<Type> v(n_reff_per_group);         // Reusable vector for transformed effects
  
  for(int i = 0; i < n_groups; ++i){
      // Reconstruct the u-vector for this group
      // Assumes u is ordered [eff1_g1..N, eff2_g1..N, ..., effk_g1..N]
      for(int j = 0; j < n_reff_per_group; ++j){
          u_group_i(j) = u(j * n_groups + i); 
      }
      
      // Solve L_cov * v = u_group_i for v
      // Use TMB's atomic solver for numerical stability
      v = atomic::solve(L_cov, u_group_i);
      
      // Add NLL components
      nll += n_reff_per_group * 0.5 * log_2pi; // k/2 * log(2pi)
      nll += logdet_L_cov;                      // log(det(L_cov))
      nll += 0.5 * v.squaredNorm();             // 0.5 * ||v||^2
  }
  
  // --- Likelihood for the Response ---
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  REPORT(betas);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  
  // Report variance components
  matrix<Type> R = L_corr * L_corr.transpose();
  
  REPORT(stdevs);
  REPORT(R); // Report the full correlation matrix
  
  // ADREPORT for standard errors
  ADREPORT(stdevs);
  ADREPORT(R);

  return nll;
}
