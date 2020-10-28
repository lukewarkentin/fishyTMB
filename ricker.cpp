#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stock);
  DATA_INTEGER(n_stocks);

  PARAMETER_VECTOR(logA); // log alpha
  PARAMETER_VECTOR(logB); // log beta
  PARAMETER_VECTOR(logSigma); // log of sigma = SD

  Type ans = 0.0;
  int N_Obs = S.size(); 
  vector<Type> LogR_Pred(N_Obs);
  vector<Type> sigma = exp(logSigma);
  vector<Type> B = exp(logB);
  
  // Ricker likelihood
  for(int i=0; i<N_Obs; ++i){
    LogR_Pred(i) = logA(stock(i)) + log(S(i)) - exp(logB(stock(i))) * S(i);
    ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stock(i)), true);
  }
  
  ADREPORT(B);
  ADREPORT(logA);
  ADREPORT(sigma);
  
  return ans;
  
}
