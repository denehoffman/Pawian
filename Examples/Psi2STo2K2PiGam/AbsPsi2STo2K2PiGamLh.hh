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
#include "Examples/Psi2STo2K2PiGam/FitParams2K2PiGam.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class Psi2STo2K2PiGamEvtList;

class AbsPsi2STo2K2PiGamLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsPsi2STo2K2PiGamLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>);
  AbsPsi2STo2K2PiGamLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>);

  /** Destructor */
  virtual ~AbsPsi2STo2K2PiGamLh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const=0;


  // Getters:
  
  double calcLogLh(const param2K2PiGam& theParamVal);
  virtual double calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const param2K2PiGam& theParamVal);

  virtual boost::shared_ptr<const Psi2STo2K2PiGamEvtList> getEventList() const {return _Psi2STo2K2PiGamEvtListPtr;}

  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal)=0;
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par)=0;
  virtual unsigned int nFitParams() =0;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal)=0;
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix)=0;


protected:

  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> _Psi2STo2K2PiGamEvtListPtr;
  FitParams2K2PiGam _fitParams2K2PiGam; 
  
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*> _evtDataVec;
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*> _evtMCVec;

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData)=0;

  virtual complex<double> chiTo2K892Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K892, double K892Mass, double K892Width);

  virtual complex<double> chiToK1K1Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK1K1, double K1Mass0, double K1Width0,  double K1Mass1, double K1Width1);

  virtual complex<double> chiTo2K_2_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_2, double K_2_Mass, double K_2_Width);

 virtual complex<double> chiToK0K2Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0K2, double K0Mass, double K0Width,  double K2Mass, double K2Width);

 virtual complex<double> chiToK0K1Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0K1, double K0Mass, double K0Width,  double K1Mass, double K1Width);

 virtual complex<double> chiToK1K2Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK1K2, double K1Mass, double K1Width,  double K2Mass, double K2Width);

  virtual  complex<double> chiToK1ToK1piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Dec, double K1Mass, double K1Width, double K1dMass, double K1dWidth);

  virtual  complex<double> chiToK0KT0KpipiKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K0dMass, double K0dWidth);

  virtual  complex<double> chiToK0KpiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0KPi, double K0Mass, double K0Width);

  virtual  complex<double> chiToKst1KpiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToKst1KPi, double Kst1Mass, double Kst1Width);

  virtual  complex<double> chiToK0KT0K1piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K1dMass, double K1dWidth);

  virtual  complex<double> chiToK1ToK0piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Dec, double K1Mass, double K1Width, double K0Mass, double K0Width);

  virtual  complex<double> chiTo2K_0_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_0, double K_0_Mass0, double K_0_Width0, double K_0_Mass1, double K_0_Width1);  

  virtual complex<double> chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif1710_k, double f980_piMass, double f980_piWidth, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif1710_k, double f980_piMass, double f980_g1, double f980_g2, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiTof980_kf0_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf0_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f0_piMass, double f0_piWidth);

  virtual complex<double> chiTof980_pif2_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif2_k, double f980_piMass, double f980_g1, double f980_g2, double f2_kMass, double f2_kWidth);

  virtual complex<double> chiTof980_kf2_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf2_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f2_piMass, double f2_piWidth);

  virtual complex<double> chiTof980f980Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980980, double f980_piMass, double f980_g1, double f980_g2);

  virtual complex<double> chiTof0_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif0_k, double f0_pi_Mass, double f0_pi_Width, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiTof0_pif2_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif2_k, double f0_pi_Mass, double f0_pi_Width, double f2_kMass, double f2_kWidth);

  virtual complex<double> chiTof2_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof2_pif0_k, double f2_pi_Mass, double f2_pi_Width, double f0_kMass, double f0_kWidth);

  virtual complex<double> chiToK_0_KToKf980KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf980K, double f980_Mass, double f980_gPiPi,  double f980_gKK, double K_0_Mass, double K_0_Width);

  virtual complex<double> chiToK_0_KToKf0KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf0K, double f0_Mass, double f0_Width, double K_0_Mass, double K_0_Width);


  virtual complex<double> chiToPi2Pi0Tof2PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double f_2_Mass, double f_2_Width);

  virtual complex<double>  chiToPi0PiToKstarKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_0_Prod, double Pi_0_Mass, double Pi_0_Width, double Kstar_Mass, double Kstar_Width);

  virtual complex<double> chiToPi2Pi0ToKstarKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double Kstar_Mass, double Kstar_Width, Spin Kstar_Spin=1);

  virtual complex<double> chiToPi0Pi0Tof980kkAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double f980_Mass, double f980_gKK, double f980_gPiPi, double Pi_0_Mass, double Pi_0_Width);

  virtual complex<double> chiToPi0Pi0Tof0Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double Pi_0_Mass, double Pi_0_Width, double f0Mass, double f0Width);


  virtual complex<double> chiToPi2Pi0Tof0PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double f0Mass, double f0Width);

  virtual complex<double> chiToK2KToK2PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K_2_Dec, double K_2prod_Mass, double K_2prod_Width, double K_2decMass, double K_2decWidth);

  virtual complex<double> chiToPi0Pi0ToK0KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double Pi_0_Mass, double Pi_0_Width, double K0Mass, double K0Width); 

  virtual  complex<double> chiToKjTof980_piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& KjProd, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& KjDec, Spin Kj, double KjMass, double KjWidth, double f980_Mass, double f980_gKK, double f980_gPiPi);

  virtual  complex<double> chiToKjTof0_piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& KjProd, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& KjDec, Spin Kj, double KjMass, double KjWidth, double f0Mass, double f0Width);

private:
  

};

#endif
