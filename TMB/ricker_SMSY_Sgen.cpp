#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  // data: 
  DATA_VECTOR(S); // spawners
  DATA_VECTOR(logR); // log(recruits)
  DATA_IVECTOR(stock); // stock ID. a vector of integers must be called in as DATA_IVECTOR(). Also MUST START AT 0 if indexing
  DATA_INTEGER(n_stocks); // number of stocks
  
  // parameters:
  PARAMETER_VECTOR(logA); // log(alpha) for each stock
  PARAMETER_VECTOR(logB); // log(beta) for each stock
  PARAMETER_VECTOR(logSigma); // log(sigma) for each stock
  PARAMETER_VECTOR(logSgen); // log(Sgen), Sgen is the spawner abundance required to reach SMSY in 1 generation
  
  // procedures (transformed parameters):
  int n = S.size(); // make integer of number data points to loop over (length of S vector)
  vector<Type> logR_Pred(n); // declare predicted logR to use in likelihood function, of length n
  vector<Type> sigma = exp(logSigma); // don't have to declare length of local variables if they are derived from parameter inputs
  vector<Type> B = exp(logB);  // don't have to declare length of local variables if they are derived from parameter inputs
  vector<Type> SMSY(n_stocks);  
  vector<Type> logSMSY(n_stocks); 
  vector<Type> Sgen = exp(logSgen);


  Type ans = 0.0; // initialize negative log likelihood at 0.0

  // Ricker likelihood
  for(int i=0; i<n; i++){ // for loop from 0 to (n-1)
    logR_Pred(i) = logA(stock(i)) + log(S(i)) - exp(logB(stock(i))) * S(i); // Ricker equation to give predicted log(recruits)
    ans -= dnorm(logR_Pred(i), logR(i),  sigma(stock(i)), true); // calculate negative log-likelihood, true is for log = TRUE. Can also use ans += -dnorm(...)
  }
  
  // Now estimate SMSY, Sgen
  SMSY = logA*(0.5-0.07*logA)/B; // Hilborn and Walters estimation of SMSY and Sgen
  logSMSY = logA + logSgen - B * Sgen;
  vector<Type> Diff = logSMSY-log(SMSY);
  ans += -sum(dnorm(Diff, 0, 1 ));
  
  
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  ADREPORT(B);
  ADREPORT(logA);
  ADREPORT(sigma);

  return ans;

}
