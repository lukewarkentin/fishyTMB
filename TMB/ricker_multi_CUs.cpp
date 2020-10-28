#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
  // data: 
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stock); // a vector of integers must be called in as DATA_IVECTOR()

  // parameters:
  PARAMETER_VECTOR(logA); // log alphas for each CU
  PARAMETER_VECTOR(logB); // log betas for each CU
  PARAMETER_VECTOR(logSigma); // log sigma for each CU
  
  // procedures (transformed parameters):
  int n = S.size(); // get number data points to loop over (length of S vector)
  
  vector<Type> logR_Pred(n); // declare predicted logR to use in likelihood function
  vector<Type> sigma = exp(logSigma); // get CU-specific sigmas
  vector<Type> B = exp(logB); // get CU-specific betas
  
  Type ans = 0.0; // initialize negative log likelihood
  
  // Ricker likelihood
  for(int i=0; i<n; ++i){ // for loop from 0 to (n-1)
    logR_Pred(i) = logA(stock(i)) + log(S(i)) - exp(logB(stock(i))) * S(i);
    ans += -dnorm(logR_Pred(i), logR(i),  sigma(stock(i)), true); // calculate negative log-likelihood, true is for log = TRUE
  }
  
  ADREPORT(B);
  ADREPORT(logA);
  ADREPORT(sigma);
  
  return ans;
  
}
