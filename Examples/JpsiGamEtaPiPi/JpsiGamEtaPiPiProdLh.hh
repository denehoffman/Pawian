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
					 double PsiToXGamMag, double PsiToXGamPhi); 

  virtual complex<double> XToPiPiEtaAmp(Spin jX, Spin lamX, EvtData* theData, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaMag, 
					 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaPhi); 

  virtual complex<double> XToAPiAmp(Spin jX, Spin lamX, Spin jA, EvtData* theData, 
				    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiMag, 
				    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiPhi,
                                    complex<double>& dynAplus, complex<double>& dynAminus);

  virtual complex<double> XToEtaFAmp(Spin jX, Spin lamX, Spin jf, EvtData* theData, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi, 
				     complex<double> dynamicFct);
  
  complex<double> A980DecFlatte(fitParams& theParamVal, const Vector4<double> &__p4);   

  
  
  bool _etaToPiPiEtaHyp;
  bool _etaToa980PiHyp;
  bool _etaToa2_1320PiHyp;
  bool _etaToSigmaEtaHyp;
  bool _etaTof0_980EtaHyp;
  bool _etaTof2_1270EtaHyp;
  bool _eta2ToPiPiEtaHyp;
  bool _eta2Toa980PiHyp;
  bool _eta2Toa2_1320PiHyp;
  bool _eta2ToSigmaEtaHyp;
  bool _eta2Tof0_980EtaHyp;
  bool _eta2Tof2_1270EtaHyp;
  bool _f1ToPiPiEtaHyp;
  bool _f1Toa980PiHyp;
  bool _f1Toa2_1320PiHyp;
  bool _f1ToSigmaEtaHyp;
  bool _f1Tof0_980EtaHyp;
  bool _f1Tof2_1270EtaHyp;
  
  bool _usePhasespace;
  
 
  

private:

  const double _massPi0;
  const double _massKplus;
  const double _massK0;
  const double _massEta;
  
  std::pair <const double, const double> _decPairPiPi;
  std::pair <const double, const double> _decPairKK;
  std::pair <const double, const double> _decPairPi0Eta;  
  bool initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   );


  

};

#endif


