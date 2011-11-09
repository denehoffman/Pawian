#ifndef _AbsPsi2SToKpKmPiGamLh_H
#define _AbsPsi2SToKpKmPiGamLh_H

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
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/FitParamsBase.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class Psi2SToKpKmPiGamEventList;

class AbsPsi2SToKpKmPiGamLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>);
  AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>);

  /** Destructor */
  virtual ~AbsPsi2SToKpKmPiGamLh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const =0;


  // Getters:
  
  virtual double calcLogLh(fitParams& theParamVal);

  virtual double calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, fitParams& theParamVal);

  virtual boost::shared_ptr<const Psi2SToKpKmPiGamEventList> getEventList() const {return _Psi2SToKpKmPiGamEvtListPtr;}

   virtual void print(std::ostream& os) const;

protected:

  boost::shared_ptr<const Psi2SToKpKmPiGamEventList> _Psi2SToKpKmPiGamEvtListPtr;

 
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtDataVec;
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtMCVec;


  virtual complex<double> a980Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiPhi, double a980Mass, double a980Width, Spin& lamChi);

  virtual complex<double> a980FlatteAmp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, fitParams& theParamVal, Spin& lamChi);


//   virtual complex<double> a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi);

   virtual complex<double> K0_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KPhi , double K0_Mass, double K0_Width, Spin& lamChi);

virtual complex<double> K1_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_1_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_1_KPhi , double K1_Mass, double K1_Width, Spin& lamChi);

  virtual complex<double> K2_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KPhi , double K2_1400Mass, double K2_1400Width, Spin& lamChi);

//   virtual complex<double> KKPi_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi, Spin& lamChi);

//   virtual complex<double> KappaPole_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappa_K, double KappaMass, double KappaWidth, Spin& lamChi);

virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData)=0;
private:


};

#endif
