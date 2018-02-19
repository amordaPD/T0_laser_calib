/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOASYMGAUSS
#define ROOASYMGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooAsymGauss : public RooAbsPdf {
public:
  RooAsymGauss() {} ; 
  RooAsymGauss(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _m0,
	      RooAbsReal& _sigma_L,
	      RooAbsReal& _sigma_R);
  RooAsymGauss(const RooAsymGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooAsymGauss(*this,newname); }
  inline virtual ~RooAsymGauss() { }
  virtual Double_t analyticalIntegral( Int_t code, const char* rangeName=0 ) const;
protected:

  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy sigma_L ;
  RooRealProxy sigma_R ;
  
  Double_t evaluate() const ;
  Double_t ApproxErf(Double_t arg) const ;
 
private:

  ClassDef(RooAsymGauss,1) // Your description goes here...
};
 
#endif
