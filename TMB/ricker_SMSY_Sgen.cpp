#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stock);
  DATA_IVECTOR(year);
  DATA_INTEGER(n_stocks);
  //DATA_VECTOR(Scales);
  
  PARAMETER_VECTOR(logA); // log alpha
  PARAMETER_VECTOR(logB); // log beta
  PARAMETER_VECTOR(logSigma);
  PARAMETER_VECTOR(logSgen);
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector<Type> LogR_Pred(N_Obs);
  vector <Type> sigma=exp(logSigma);
  vector <Type> SMSY(n_stocks);  
  vector <Type> logSMSY(n_stocks);
  vector <Type> Sgen = exp(logSgen);
  vector <Type> B = exp(logB);
  
  // Ricker likelihood
  for(int i=0; i<N_Obs; i++){
    LogR_Pred(i) = logA(stock(i)) + log(S(i)) - exp(logB(stock(i))) * S(i);
    ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stock(i)), true);
  }
  
  // Now estimate SMSY, Sgen
  SMSY = logA*(0.5-0.07*logA)/B; // Hilborn and Walters estimation of SMSY and Sgen
  logSMSY = logA + logSgen - B * Sgen;
  vector <Type> Diff = logSMSY-log(SMSY);
  ans += -sum(dnorm(Diff, 0, 1 ));
  
  
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  
  return ans;
  
}
