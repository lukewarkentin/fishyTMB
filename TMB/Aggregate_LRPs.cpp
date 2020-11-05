#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  // data 
  DATA_VECTOR(S); // spawners
  DATA_VECTOR(logR); // log(recruits)
  DATA_IVECTOR(stock); // stock / CU name
  DATA_IVECTOR(yr); // year
  DATA_INTEGER(n_stocks); // number of stocks / CUs
  DATA_INTEGER(Mod_Yr_0); // first year of SR data
  DATA_INTEGER(Mod_Yr_n); // last year of SR data

  // parameters
  PARAMETER_VECTOR(logA); 
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);
  PARAMETER_VECTOR(logSgen);
  PARAMETER(B_0); // binomial distribution parameter
  PARAMETER(B_1); // binomial distribution parameter
  
  // procedures (transformed parameters):
  int n = S.size();  // make integer of number data points to loop over (length of S vector)
  vector<Type> logR_Pred(n); // vector of predicted log(recruits)
  vector <Type> sigma=exp(logSigma); 
  vector <Type> B = exp(logB);
  vector <Type> SMSY(n_stocks);  // stock-specific SMSY values
  vector <Type> logSMSY(n_stocks); // stock-specific log(SMSY) values
  vector <Type> Sgen = exp(logSgen); 
  
  Type ans=0.0; // initialize log-likelihood at 0.0
  
  // Ricker likelihood
  for(int i=0; i<n; i++){
    logR_Pred(i) = logA(stock(i)) + log(S(i)) - exp(logB(stock(i))) * S(i); // Ricker equation to give predicted log(recruits)
    ans += -dnorm(logR_Pred(i), logR(i),  sigma(stock(i)), true); // calculate negative log-likelihood, true is for log = TRUE.
  }
  
  // Now estimate SMSY, Sgen
  SMSY = logA*(0.5-0.07*logA)/B;
  logSMSY = logA + logSgen - B * Sgen;
  vector <Type> Diff = logSMSY-log(SMSY);
  ans += -sum(dnorm(Diff, 0, 1 ));
  
  // go through ETS for each year and see how many stocks (what is ETS?)
  // are above their benchmark
  int Logistic_Mod_Yrs = Mod_Yr_n - Mod_Yr_0 + 1;
  vector <Type> N_Above_LRP(Logistic_Mod_Yrs);
  N_Above_LRP.setZero();
  vector <Type> Agg_Abund(Logistic_Mod_Yrs);
  Agg_Abund.setZero();
  
  for(int i=0; i<n; ++i){
    if(yr(i) >= Mod_Yr_0 && yr(i) <= Mod_Yr_n){
      //check if ETS above LRP
      if(S(i) > Sgen(stock(i))){
        N_Above_LRP(yr(i)-Mod_Yr_0) += 1;
      }
      //add to aggregate abund
      Agg_Abund(yr(i)-Mod_Yr_0) += S(i);
    } // end if model year
  } // end for loop over obs
  
  // Now add logistic model likelihood 
  vector<Type> LogitP(Logistic_Mod_Yrs);
  vector<Type> N(Logistic_Mod_Yrs);
  // Scale down Agg Abund
  //Type Agg_Sum = sum(Agg_Abund);
  //Type Agg_Mean = Agg_Sum/N_Obs;
  
  for(int i=0; i<Logistic_Mod_Yrs; i++){
    LogitP(i) = B_0 + B_1*Agg_Abund(i);
    N(i) = n_stocks;
  }
  ans += -sum(dbinom_robust(N_Above_LRP, N, LogitP, true));
  
  //Get final BM
  //Type Agg_BM = ((log(0.95 / 0.05) - B_0)*Agg_Mean)/(B_1);
  Type Agg_BM = (log(0.95 / 0.05) - B_0)/(B_1);
  
  ADREPORT(B);
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  REPORT(N_Above_LRP);
  REPORT(Agg_Abund);
  ADREPORT(Agg_BM);
  
  
  return ans;
  
}
