#ifndef _JpsiGamEtaPiPiProdLhNew_H
#define _JpsiGamEtaPiPiProdLhNew_H

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

#include "PwaUtils/AbsLhNew.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"

#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;

class JpsiGamEtaPiPiProdLhNew : public AbsLhNew{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiProdLhNew(boost::shared_ptr<const EvtDataBaseList>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);
  JpsiGamEtaPiPiProdLhNew(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);

  /** Destructor */
  virtual ~JpsiGamEtaPiPiProdLhNew();

  virtual AbsLhNew* clone_() const {
    return new JpsiGamEtaPiPiProdLhNew(_evtListPtr, _hypVec, _theStatesPtr);
  }
  
  virtual double calcEvtIntensity( EvtData* theData, fitParamsNew& theParamVal);
  
  //Getters:
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  
  
protected:
  const std::vector<std::string> _hypVec;
  std::vector<std::string> _GammaEtaHyps;
  std::vector<std::string> _GammaEta2Hyps; 
  std::vector<std::string> _GammaF1Hyps;
  boost::shared_ptr<JpsiGamEtaPiPiStates> _theStatesPtr; 
  bool _usePhasespace;  
 
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _etaDecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _f1DecAmpMap;

  virtual complex<double> psiToXGammaAmp(Spin Minit, Spin jX, Spin lamX, Spin lamGamma, EvtData* theData, 
                                         double PsiToXGamMag, double PsiToXGamPhi);
private:

  const std::string _GammaEtaKey;
  const std::string _GammaEta2Key;
  const std::string _GammaF1Key;
  const std::string _EtaKey;
  const std::string _Eta2Key;
  const std::string _F1Key;

  const double _massPi0;
  const double _massKplus;
  const double _massK0;
  const double _massEta;
  
  std::pair <const double, const double> _decPairPiPi;
  std::pair <const double, const double> _decPairKK;
  std::pair <const double, const double> _decPairPi0Eta;  

  void initializeHypothesis();
  void getDefaultLamLamParams(std::vector<std::string>& hyps, std::vector< boost::shared_ptr<const JPClamlam> > lamLamAmps, fitParamsNew& fitVal, fitParamsNew& fitErr);

  

};

#endif


