/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "RooFit.h"
#include <math.h>

#include "Riostream.h" 

#include "RooExpGauss.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooRealVar.h"
#include "RooMath.h"
#include <math.h> 
#include "TMath.h" 

ClassImp(RooExpGauss) 

 RooExpGauss::RooExpGauss(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _sigma,
                        RooAbsReal& _alpha) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   m0("m0","m0",this,_m0),
   sigma("sigma","sigma",this,_sigma),
   alpha("alpha","alpha",this,_alpha)
 { 
 } 


 RooExpGauss::RooExpGauss(const RooExpGauss& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   m0("m0",this,other.m0),
   sigma("sigma",this,other.sigma),
   alpha("alpha",this,other.alpha)
 { 
 } 



 Double_t RooExpGauss::evaluate() const 
 {
   
  
  Double_t t = (m-m0)/sigma;
  if (alpha < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)alpha);

  if (t >= -absAlpha) {
    return exp(-0.5*t*t);
  }
  else {
    Double_t a = exp(0.5*absAlpha*absAlpha+absAlpha*t);

    return a;
  }
 } 



Double_t RooExpGauss::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static const double sqrtPiOver2 = 1.2533141373;
  static const double sqrt2 = 1.4142135624;

  R__ASSERT(code==1);
  double result = 0.0;

  
  double sig = fabs((Double_t)sigma);
  
  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;
  
  if(alpha < 0) {
    double tmp = tmin;
    tmin = -tmax;
    tmax = -tmp;
  }

  double absAlpha = fabs((Double_t)alpha);
  
  if( tmin >= -absAlpha ) {
    result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                - ApproxErf(tmin/sqrt2) );
  } //� tutta gaussian nell'intervallo dato
  else if( tmax <= -absAlpha ) { //� tutta esponenziale nell'intervallo dato    
    result += sig/absAlpha*exp(0.5*absAlpha*absAlpha)*(exp(absAlpha*tmax)-exp(absAlpha*tmin))  ;
 
  }
  else { //� un p� gaussian, un p� esponenziale
 
    
    double term1 = 0.0;
    term1 = sig/absAlpha*exp(0.5*absAlpha*absAlpha)*(exp(absAlpha*absAlpha)-exp(absAlpha*tmin))  ;
    double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                     - ApproxErf(-absAlpha/sqrt2) );
    
    
    result += term1 + term2;
  }
  
  return result;
}

Double_t RooExpGauss::ApproxErf(Double_t arg) const 
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return RooMath::erf(arg);
}
