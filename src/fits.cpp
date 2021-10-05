// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data input
DATA_VECTOR(w)
DATA_VECTOR(fecundity)
DATA_INTEGER(nobs)

PARAMETER(lalpha)
PARAMETER(beta)

vector<Type>fits(nobs);

Type alpha = exp(lalpha);


for(int i=0;i<nobs;i++){ // Start time loop
 fits(i) = alpha*pow(w(i), beta);
}


// Likelihood
Type ans = 0.0;

for(int i=0;i<nobs;i++){

ans += -dnorm(fecundity(i), fits(i), Type(1.0), TRUE);

}
//


REPORT(fits)


ADREPORT(fits)



  return ans;
}
