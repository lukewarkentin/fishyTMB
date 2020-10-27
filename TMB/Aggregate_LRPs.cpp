#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);
  DATA_INTEGER(N_Stks);
  DATA_INTEGER(Mod_Yr_0);
  DATA_INTEGER(Mod_Yr_n);
  //DATA_VECTOR(Scales);
  
  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);
  PARAMETER_VECTOR(logSgen);
  PARAMETER(B_0);
  PARAMETER(B_1);
  
  Type ans=0.0;
  int N_Obs = S.size(); 
  vector<Type> LogR_Pred(N_Obs);
  vector <Type> sigma=exp(logSigma);
  vector <Type> SMSY(N_Stks);  
  vector <Type> LogSMSY(N_Stks);
  vector <Type> Sgen = exp(logSgen);
  vector <Type> B = exp(logB);
  
  // Ricker likelihood
  for(int i=0; i<N_Obs; i++){
    LogR_Pred(i) = logA(stk(i)) + log(S(i)) - exp(logB(stk(i))) * S(i);
    ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stk(i)), true);
  }
  
  // Now estimate SMSY, Sgen
  SMSY = logA*(0.5-0.07*logA)/B;
  LogSMSY = logA + logSgen - B * Sgen;
  vector <Type> Diff = LogSMSY-log(SMSY);
  ans += -sum(dnorm(Diff, 0, 1 ));
  
  // go through ets for each year and see how many stocks
  // are above their benchmark
  int Logistic_Mod_Yrs = Mod_Yr_n - Mod_Yr_0 + 1;
  vector <Type> N_Above_LRP(Logistic_Mod_Yrs);
  N_Above_LRP.setZero();
  vector <Type> Agg_Abund(Logistic_Mod_Yrs);
  Agg_Abund.setZero();
  
  for(int i=0; i<N_Obs; ++i){
    if(yr(i) >= Mod_Yr_0 && yr(i) <= Mod_Yr_n){
      //check if ETS above LRP
      if(S(i) > Sgen(stk(i))){
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
  
  for(int i=0; i<Logistic_Mod_Yrs; ++i){
    LogitP(i) = B_0 + B_1*Agg_Abund(i);
    N(i) = N_Stks;
  }
  ans += -sum(dbinom_robust(N_Above_LRP, N, LogitP, true));
  
  //Get final BM
  //Type Agg_BM = ((log(0.95 / 0.05) - B_0)*Agg_Mean)/(B_1);
  Type Agg_BM = (log(0.95 / 0.05) - B_0)/(B_1);
  
  
  ADREPORT(SMSY);
  ADREPORT(Sgen);
  REPORT(N_Above_LRP);
  REPORT(Agg_Abund);
  ADREPORT(Agg_BM);
  
  
  return ans;
  
}
