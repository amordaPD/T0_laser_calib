/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "RooFit.h"
#include <math.h>

#include "Riostream.h" 

#include "RooAsymGauss.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooRealVar.h"
#include "RooMath.h"
#include <math.h> 
#include "TMath.h" 

ClassImp(RooAsymGauss) 

 RooAsymGauss::RooAsymGauss(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _sigma_L,
                        RooAbsReal& _sigma_R) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   m0("m0","m0",this,_m0),
   sigma_L("sigma_L","sigma_L",this,_sigma_L),
   sigma_R("sigma_R","sigma_R",this,_sigma_R)
 { 
 } 


 RooAsymGauss::RooAsymGauss(const RooAsymGauss& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   m0("m0",this,other.m0),
   sigma_L("sigma_L",this,other.sigma_L),
   sigma_R("sigma_R",this,other.sigma_R)
 { 
 } 



 Double_t RooAsymGauss::evaluate() const 
 {
   
   static const double piG = 3.1415926535 ;
   double sig_L = fabs((Double_t)sigma_L);
   double sig_R = fabs((Double_t)sigma_R);
   double norm = sqrt(2.0)/(sqrt(piG)*(sig_L+sig_R));
   Double_t t = (m-m0);
   

  if (t >= 0) {
    return norm*exp(-0.5*t*t/(sigma_R*sigma_R));
  }
  else {
    return norm*exp(-0.5*t*t/(sigma_L*sigma_L));
  }
 } 



Double_t RooAsymGauss::analyticalIntegral(Int_t code, const char* rangeName) const
{
  static const double sqrtPiOver2 = 1.2533141373;
  static const double sqrt2 = 1.4142135624;
  static const double piG = 3.1415926535 ;
  R__ASSERT(code==1);
  double result = 0.0;

  
  double sig_L = fabs((Double_t)sigma_L);
  double sig_R = fabs((Double_t)sigma_R);
  double norm = sqrt(2.0)/(sqrt(piG)*(sig_L+sig_R));
  
  double tmin = (m.min(rangeName)-m0);
  double tmax = (m.max(rangeName)-m0);
  
  
  if( tmin >= 0 ) {
    result += norm*sqrt(piG)*0.5*(   ApproxErf(tmax/(sig_R*sqrt2))
				    - ApproxErf(tmin/(sig_R*sqrt2)) );
  } 
  else if( tmax < 0 ) { 
    result += norm*sqrt(piG)*0.5*(   ApproxErf(tmax/(sig_L*sqrt2))
				    - ApproxErf(tmin/(sig_L*sqrt2)) );
 
  }
  else { 
 
    
    double term1 = 0.0;
    term1 = norm*sqrt(piG)*0.5*(- ApproxErf(tmin/(sig_L*sqrt2)) );
    
    double term2 = 0.0;
    term2 = norm*sqrt(piG)*0.5*(ApproxErf(tmax/(sig_R*sqrt2)));
        
    result += term1 + term2;
  }
  
  return result;
}

Double_t RooAsymGauss::ApproxErf(Double_t arg) const 
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return RooMath::erf(arg);
}
