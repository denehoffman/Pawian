#ifndef _JpsiGamEtaPiPiProdLh_H
#define _JpsiGamEtaPiPiProdLh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"

#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"


class JpsiGamEtaPiPiProdLh : public AbsLh{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiProdLh(boost::shared_ptr<const EvtDataBaseList>, const std::map<const std::string, bool>& hypMap);
  JpsiGamEtaPiPiProdLh(boost::shared_ptr<AbsLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~JpsiGamEtaPiPiProdLh();

  virtual AbsLh* clone_() const {
    return new JpsiGamEtaPiPiProdLh(_evtListPtr, _hypMap);
  }
  
  virtual double calcEvtIntensity( EvtData* theData, fitParams& theParamVal);
  
  //Getters:
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;
  
  
protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams&  theParamVal, EvtData* theData);  

  virtual complex<double> psiToXGammaAmp(Spin Minit, Spin jX, Spin lamX, Spin lamGamma, EvtData* theData, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToXGamMag, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToXGamPhi); 

  virtual complex<double> XToPiPiEtaAmp(Spin jX, Spin lamX, EvtData* theData, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaMag, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaPhi); 

  virtual complex<double> XToAPiAmp(Spin jX, Spin lamX, Spin jA, EvtData* theData, 
				    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiMag, 
				    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiPhi,
                                    complex<double>& dynAplus, complex<double>& dynAminus);

  complex<double> A980DecFlatte(fitParams& theParamVal, const Vector4<double> &__p4);   

  virtual void calcEtaGammaAmp( EvtData* theData, 
				std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaGamMag, 
				std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaGamPhi, 
				complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm);
  

  
  
  bool _etaToPiPiEtaHyp;
  bool _etaToa980PiHyp;
  bool _etaToa2_1320PiHyp;
  bool _eta2ToPiPiEtaHyp;
  bool _eta2Toa980PiHyp;
  bool _eta2Toa2_1320PiHyp;
  bool _f1ToPiPiEtaHyp;
  bool _f1Toa980PiHyp;
  
  bool _usePhasespace;
  
  std::map<const std::string, bool> _hypMap;

 
  

private:
  
  bool initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   );


  

};

#endif


