#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLhNew.hh"
#include "PwaUtils/PsiToXGamStates.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;

class PsiProdBaseLhNew : public AbsLhNew {

public:
  PsiProdBaseLhNew(boost::shared_ptr<const EvtDataBaseListNew>, 
		   const std::vector<std::string>& hypVec, 
		   boost::shared_ptr<PsiToXGamStates> theStates);
  PsiProdBaseLhNew(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, 
		   boost::shared_ptr<PsiToXGamStates> theStates);
  virtual ~PsiProdBaseLhNew();

  virtual AbsLhNew* clone_() const = 0;
  virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  
protected:
  const std::vector<std::string> _hypVec;
  std::vector<std::string> _GammaEtaHyps;
  std::vector<std::string> _GammaEta2Hyps; 
  std::vector<std::string> _GammaF0Hyps;
  std::vector<std::string> _GammaF1Hyps;
  std::vector<std::string> _GammaF2Hyps;
  std::map< std::string,std::vector<std::string> > _hypMap; 
  boost::shared_ptr<PsiToXGamStates> _theStatesPtr;
 
  bool _usePhasespace;  
  const std::string _GammaEtaKey;
  const std::string _GammaEta2Key;
  const std::string _GammaF0Key;
  const std::string _GammaF1Key;
  const std::string _GammaF2Key;
  const std::string _EtaKey;
  const std::string _Eta2Key;
  const std::string _F0Key;
  const std::string _F1Key;
  const std::string _F2Key;
  const std::string _phasespaceKey;
 
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _etaDecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _f0DecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _f1DecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _f2DecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _eta2DecAmpMap;
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >  _allAmpMap; 


  virtual complex<double> psiToXGammaAmp(Spin Minit, Spin lamX, Spin lamGamma,
					 EvtDataNew* theData, 
                                         double PsiToXGamMag, double PsiToXGamPhi);
  
  virtual void checkParamVariation(fitParamsNew& theParamVal);
 
private:
  void getDefaultLamLamParams(const std::vector<std::string>& hyps, 
			      std::vector< boost::shared_ptr<const JPClamlam> > lamLamAmps, 
			      fitParamsNew& fitVal, fitParamsNew& fitErr);
  void initializeHypothesis();
};
