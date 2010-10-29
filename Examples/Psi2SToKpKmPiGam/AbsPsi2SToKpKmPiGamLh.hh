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
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class Psi2SToKpKmPiGamEventList;
class Psi2SToKpKmPiGamStates;

class AbsPsi2SToKpKmPiGamLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>, boost::shared_ptr<const Psi2SToKpKmPiGamStates>);
  AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>);

  /** Destructor */
  virtual ~AbsPsi2SToKpKmPiGamLh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const =0;


  // Getters:
  
  double calcLogLh(const Psi2SToKpKmPiGamData::fitParamVal& theParamVal);
  virtual double calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal);

  virtual boost::shared_ptr<const Psi2SToKpKmPiGamEventList> getEventList() const {return _Psi2SToKpKmPiGamEvtListPtr;}
  virtual boost::shared_ptr<const Psi2SToKpKmPiGamStates> getPsi2SToKpKmPiGamStates() const {return _Psi2SToKpKmPiGamStatesPtr;}

  virtual void setMnUsrParams(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal)=0;

  virtual int setFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const=0;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(Psi2SToKpKmPiGamData::fitParamVal& theParamVal) const;

protected:

  boost::shared_ptr<const Psi2SToKpKmPiGamEventList> _Psi2SToKpKmPiGamEvtListPtr;
  boost::shared_ptr<const Psi2SToKpKmPiGamStates> _Psi2SToKpKmPiGamStatesPtr;
  
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtDataVec;
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*> _evtMCVec;

  virtual int setFitParamValDec(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const; 

  virtual int setFitParamValMass(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const; 

  virtual int setFitParamFlattea980Mass(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const; 

  virtual complex<double> a980Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980Mass, double a980Width, Spin& lamChi);

  virtual complex<double> a980FlatteAmp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980FlatteMass, double a980FlatgKK, double a980FlatgEtaPi, Spin& lamChi);


  virtual complex<double> a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi);

  virtual complex<double> K892Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K, double K892Mass, double K892Width, Spin& lamChi);

  virtual complex<double> K0_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K, double K0_1400Mass, double K0_1400Width, Spin& lamChi);

  virtual complex<double> K1_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_1_K, double K1_1400Mass, double K1_1400Width, Spin& lamChi);

  virtual complex<double> K2_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_2_K, double K2_1400Mass, double K2_1400Width, Spin& lamChi);

  virtual void setMnUsrParamsDec(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key);

  virtual void setMnUsrParamsMass(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key);

  virtual void setMnUsrParamsFlattea980Mass(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key);

  virtual void checkFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& fitVal);

private:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData)=0;

};

#endif
