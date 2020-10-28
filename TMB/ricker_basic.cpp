#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
  // data: 
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  
  // parameters:
  PARAMETER(logA); // log alpha
  PARAMETER(logB); // log beta
  PARAMETER(logSigma); // log of sigma = SD

  // procedures (transformed parameters):
  int n = S.size(); // get number data points to loop over (length of S vector)
  
  vector<Type> logR_Pred(n); // declare predicted logR to use in likelihood function
  Type sigma = exp(logSigma); 
  Type B = exp(logB);
  
  Type ans = 0.0; // initialize negative log likelihood
  
  // Ricker likelihood
  for(int i=0; i<n; ++i){ // for loop from 0 to (n-1)
    logR_Pred(i) = logA + log(S(i)) - exp(logB) * S(i);
    ans += -dnorm(logR_Pred(i), logR(i),  sigma, true); // calculate negative log-likelihood, true is for log = TRUE
  }
  
  ADREPORT(B);
  ADREPORT(logA);
  ADREPORT(sigma);
  
  return ans;
  
}
