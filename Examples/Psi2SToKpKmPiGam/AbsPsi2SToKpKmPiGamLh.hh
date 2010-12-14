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
#include "Examples/Psi2SToKpKmPiGam/FitParamsKpKmPiGam.hh"
#include "PwaUtils/DataUtils.hh"

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
  
  double calcLogLh(const paramKpKmPiGam& theParamVal);
  virtual double calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const paramKpKmPiGam& theParamVal);

  virtual boost::shared_ptr<const Psi2SToKpKmPiGamEventList> getEventList() const {return _Psi2SToKpKmPiGamEvtListPtr;}

  virtual void setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal)=0;

  virtual int setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par) =0;
  virtual unsigned int nFitParams() =0;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(paramKpKmPiGam& theParamVal)=0;
  virtual void dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix)=0;
protected:

  boost::shared_ptr<const Psi2SToKpKmPiGamEventList> _Psi2SToKpKmPiGamEvtListPtr;
  FitParamsKpKmPiGam _fitParamsKpKmPiGam; 
 
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtDataVec;
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtMCVec;


  virtual complex<double> a980Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980Mass, double a980Width, Spin& lamChi);

  virtual complex<double> a980FlatteAmp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980FlatteMass, double a980FlatgKK, double a980FlatgEtaPi, Spin& lamChi);


  virtual complex<double> a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi);

  virtual complex<double> K892Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K, double K892Mass, double K892Width, Spin& lamChi);

  virtual complex<double> K0_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K, double K0_1400Mass, double K0_1400Width, Spin& lamChi);

  virtual complex<double> K1_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_1_K, double K1_1400Mass, double K1_1400Width, Spin& lamChi);

  virtual complex<double> K2_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_2_K, double K2_1400Mass, double K2_1400Width, Spin& lamChi);

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData)=0;
private:


};

#endif
