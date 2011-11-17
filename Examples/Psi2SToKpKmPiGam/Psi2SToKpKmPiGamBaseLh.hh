#ifndef _Psi2SToKpKmPiGamBaseLh_H
#define _Psi2SToKpKmPiGamBaseLh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEvtList.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

using namespace ROOT::Minuit2;

class Psi2SToKpKmPiGamBaseLh : public AbsLh{

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<const EvtDataBaseList>, const std::map<const std::string, bool>& hypMap);
  Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<AbsLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Psi2SToKpKmPiGamBaseLh();

  virtual AbsLh* clone_() const {
    return new Psi2SToKpKmPiGamBaseLh(_evtListPtr, _hypMap);
  }

  virtual double calcEvtIntensity(EvtData* theData, fitParams& theParamVal);
  // Getters:

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;

protected:
  bool _K0_1430Hyp;
  bool _K1_1410Hyp;
  bool _K2_1430Hyp;
  bool _K1_1680Hyp;
  bool _KKPi_Hyp;
  std::map<const std::string, bool> _hypMap;
  virtual  complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData);

  virtual  complex<double> calcDecAmp(complex<double>& inAmp, Spin lamChi, fitParams& theParamVal, EvtData* theData);

  virtual complex<double> a980Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiPhi, double a980Mass, double a980Width, Spin& lamChi);

  virtual complex<double> a980FlatteAmp(EvtData* theData, fitParams& theParamVal, Spin& lamChi);


//   virtual complex<double> a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi);

   virtual complex<double> K0_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KPhi , double K0_Mass, double K0_Width, Spin& lamChi);

virtual complex<double> K1_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_1_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_1_KPhi , double K1_Mass, double K1_Width, Spin& lamChi);

  virtual complex<double> K2_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KPhi , double K2_1400Mass, double K2_1400Width, Spin& lamChi);

//   virtual complex<double> KKPi_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi, Spin& lamChi);

//   virtual complex<double> KappaPole_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappa_K, double KappaMass, double KappaWidth, Spin& lamChi);


private:


};

#endif
