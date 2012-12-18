#ifndef _XDecAmpBase_H
#define _XDecAmpBase_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiStates.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

class XDecAmpBase : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  XDecAmpBase(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates, Spin spinX);

  /** Destructor */
  virtual ~XDecAmpBase();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParamsNew& theParamVal);

  virtual void updateFitParams(fitParamsNew& theParamVal);

protected:
  const std::string _piPiEtaKey;
  bool _piPiEtaHyp;
  const std::string _a980piKey;
  bool _a980piHyp;
  const std::string _sigmaEtaKey;
  bool _sigmaEtaHyp;
  const std::string _f980etaKey;
  bool _f980etaHyp;
  const std::string _a2_1320piKey;
  bool _a2_1320piHyp;
  const std::string _f2_1270etaKey;
  bool _f2_1270etaHyp;
  const std::string _pipiSetaKey;
  bool _pipiSetaHyp;
  const std::string _xBWKey;
  bool _massIndependent;

  const double _massPi;
  const double _massPi0;
  const double _massKplus;
  const double _massK0;
  const double _massEta;
  std::pair <const double, const double> _decPairKK;
  std::pair <const double, const double> _decPairPiEta;
  std::pair <const double, const double> _decPairPi0Eta;
  std::pair <const double, const double> _decPairPiPi;
  std::pair <const double, const double> _decPairPi0Pi0;
  boost::shared_ptr<JpsiGamEtaPiPiStates> _theStatesPtr;
  bool _recalculatef2_1270;
  bool _recalculatea2_1320;
  bool _recalculatepipiSeta;

  complex<double> XToPiPiEtaAmp(Spin lamX, EvtDataNew* theData,
				std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaMag,
				std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaPhi);
  complex<double> XToAPiFlatteAmp(Spin lamX, EvtDataNew* theData,
				  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToA980PiMag,
				  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToA980PiPhi,
				  double a0_980Mass, double a0_980gPiEta,double a0_980gKK);
  complex<double> XToFEtaFlatteAmp(Spin lamX, EvtDataNew* theData,
				   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToFEtaFlatteMag,
				   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToFEtaFlattePhi,
				   double f0_980Mass, double f0_980gPiPi,double f0_980gKK);
  complex<double> XToEtaFAmp(Spin lamX, Spin jf, EvtDataNew* theData, 
			     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
			     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi,
			     double fMass, double fWidth);
  complex<double> XToPiPiSEtaAmp(Spin lamX, EvtDataNew* theData);
  complex<double> XToAPiBWAmp(Spin lamX, Spin jA, EvtDataNew* theData, 
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiMag,
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiPhi, 
			      double aMass, double aWidth);

  virtual void initialize();

private:
  boost::shared_ptr<FVectorPiPiS> _pipiSFVec;  
  std::map<int, std::map<Spin, complex<double> > > _cachedAmpa2_1320;
  std::map<int, std::map<Spin, complex<double> > > _cachedAmpf2_1270;
  std::map<int, std::map<Spin, complex<double> > > _cachedAmppipiSeta;

  std::vector<std::string> _enabledAmpKeys;
  std::vector<std::string> _enabledMassKeys;
  std::vector<std::string> _enabledWidthKeys;
  std::vector<std::string> _enabledFactorKeys;
  std::map<std::string, double> _currentMassMap;
  std::map<std::string, double> _currentWidthMap;
  std::map<std::string, double> _currentgFactorMap;
  std::map<std::string, double> _currentbFactorMap;
  std::map<std::string, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > _currentParamMagMap;
  std::map<std::string, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > _currentParamPhiMap;
};

#endif


