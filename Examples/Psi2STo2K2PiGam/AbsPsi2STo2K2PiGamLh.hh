#ifndef _AbsPsi2STo2K2PiGamLh_H
#define _AbsPsi2STo2K2PiGamLh_H

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
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class Psi2STo2K2PiGamEvtList;

class AbsPsi2STo2K2PiGamLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsPsi2STo2K2PiGamLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, boost::shared_ptr<const Psi2STo2K2PiGamStates>);
  AbsPsi2STo2K2PiGamLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>);

  /** Destructor */
  virtual ~AbsPsi2STo2K2PiGamLh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const=0;


  // Getters:
  
  double calcLogLh(const Psi2STo2K2PiGamData::fitParamVal& theParamVal);
  virtual double calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const Psi2STo2K2PiGamData::fitParamVal& theParamVal);

  virtual boost::shared_ptr<const Psi2STo2K2PiGamEvtList> getEventList() const {return _Psi2STo2K2PiGamEvtListPtr;}
  virtual boost::shared_ptr<const Psi2STo2K2PiGamStates> getPsi2STo2K2PiGamStates() const {return _Psi2STo2K2PiGamStatesPtr;}

  virtual void setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal);

  virtual int setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const;
  virtual void dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const;


protected:

  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> _Psi2STo2K2PiGamEvtListPtr;
  boost::shared_ptr<const Psi2STo2K2PiGamStates> _Psi2STo2K2PiGamStatesPtr;
  
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*> _evtDataVec;
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*> _evtMCVec;

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

  virtual complex<double> chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData)=0;

  virtual int setFitParamValDec(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const; 

  virtual int setFitParamValMass(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const; 
  
  virtual int setFitParamFlatteMass(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const;
  
    virtual void setMnUsrParamsDec(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key);

  virtual void setMnUsrParamsMass(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key);

  virtual void setMnUsrParamsFlatteMass(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key);

  virtual void checkFitParamVal(Psi2STo2K2PiGamData::fitParamVal& fitVal);

  virtual complex<double> chiTo2K892Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K892, double K892Mass, double K892Width);

  virtual complex<double> chiToK1K1Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK1K1, double K1Mass0, double K1Width0,  double K1Mass1, double K1Width1);

  virtual complex<double> chiTo2K_2_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_2, double K_2_Mass, double K_2_Width);

  virtual  complex<double> chiToK1ToK1piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1ToK1dPi, double K1Mass, double K1Width, double K1dMass, double K1dWidth);


  virtual  complex<double> chiToK0KT0KpipiKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K0dMass, double K0dWidth);

  virtual  complex<double> chiToK0KT0K1piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K1dMass, double K1dWidth);

  virtual  complex<double> chiToK1ToK0piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1ToK0Pi, double K1Mass, double K1Width, double K0Mass, double K0Width);

  virtual  complex<double> chiTo2K_0_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_0, double K_0_Mass0, double K_0_Width0, double K_0_Mass1, double K_0_Width1);  

  virtual complex<double> chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif1710_k, double f980_piMass, double f980_piWidth, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif1710_k, double f980_piMass, double f980_g1, double f980_g2, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiTof980_kf0_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf0_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f0_piMass, double f0_piWidth);

  virtual complex<double> chiTof980f980Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980980, double f980_piMass, double f980_g1, double f980_g2);

  virtual complex<double> chiTof0_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif0_k, double f0_pi_Mass, double f0_pi_Width, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiToK_0_KToKf980KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf980K, double f980_Mass, double f980_gPiPi,  double f980_gKK, double K_0_Mass, double K_0_Width);


private:
  

};

#endif
