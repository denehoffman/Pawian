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
//#include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiDynamicModelParams.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"


class FitParamErrorMatrix;

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
  
  virtual double calcEvtIntensity(EvtData* theData, fitParams& theParamVal);
  
  //Getters:
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;
  
  void useCommonProductionPhase(bool commonPhase){_useCommonProductionPhase=commonPhase;}
  void massIndependentFit(bool massIndep){ _massIndependentFit=massIndep;  }
  
  bool calcComponentIntensity(EvtData* theData, fitParams& theParamVal,  FitParamErrorMatrix& theErrMatrix, std::string component, std::pair<double, double> &intensity);
  void dumpComponentIntensity(std::ostream &os, fitParams& theParams, FitParamErrorMatrix& theErrMatrix);
  

  
protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams&  theParamVal, EvtData* theData);  
  
  virtual complex<double> etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, EvtData* theData, 
				      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaMag, 
				      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaPhi);
  
  virtual complex<double> f0GammaAmp(Spin Minit, Spin Mgamma, EvtData* theData, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdMag, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdPhi,
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag,  
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi,
				     bool useCommonProductionPhase, dynamicModelParams& dynModPars);
  /*
  virtual complex<double> f2GammaAmp(Spin Minit, Spin Mgamma, Spin fJSpin, EvtData* theData, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdMag, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdPhi,
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecMag,  
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecPhi,
				     double mass, double width, bool useCommonProductionPhase, dynamicModelParams& dynModPars);
  virtual complex<double> etaToPhiPhiTo4KAmp(EvtData* theData);
  virtual complex<double> f0ToPhiPhiTo4KAmp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag,
		  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi, dynamicModelParams& dynModPars);
  complex<double> f2ToPhiPhiTo4KAmp(EvtData* theData, Spin f2Lambda,
		  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecMag ,
		  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecPhi,
		  dynamicModelParams& dynModPars);
  complex<double> phiphiTo4KAmp(EvtData* theData, Spin lambdaPhi1, Spin lambdaPhi2);
   */

  bool  initializeFlatteModel(dynamicModelParams &theDynModel, const Vector4<double> &fv2Phi, double mass, double gPhiPhi, double gKK);
  bool  initializeBreitWignerModel(dynamicModelParams &theDynModel,const Vector4<double> &fv2Phi, double mass, double width, dynamicModelParams::enumDynamicModel theModel);
  
  void calcF11285Amp(EvtData* theData, fitParams& theParamVal, complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm,
		  dynamicModelParams::enumDynamicModel theModel);
  void calcEta1295Amp(EvtData* theData, fitParams& theParamVal, complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm,
		  dynamicModelParams::enumDynamicModel theModel);
  void calcEta1405Amp(EvtData* theData, fitParams& theParamVal, complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm,
		  dynamicModelParams::enumDynamicModel theModel);

  /*
  void calcF22300Amp(EvtData* theData, fitParams& theParamVal,complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
		  dynamicModelParams::enumDynamicModel theModel);
  void calcE21870Amp(EvtData* theData, fitParams& theParamVal, complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
		  dynamicModelParams::enumDynamicModel theModel);
  void calcF1Amp(EvtData* theData, fitParams& theParamVal, complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
		  dynamicModelParams::enumDynamicModel theModel);
  */
  
  
  bool _f11285Hyp;
  bool _eta1295Hyp;
  bool _eta1405Hyp;
  
  bool _usePhasespace;
  bool _useCommonProductionPhase;
  bool _massIndependentFit;
  
  std::map<const std::string, bool> _hypMap;

 
  

private:
  
  void calcAmpError(EvtData* theData, fitParams& theParamVal, FitParamErrorMatrix& theErrMatrix,
 		       dynamicModelParams::enumDynamicModel theModel, 
		     boost::function<void(EvtData* , fitParams&, complex<double>&, 
					  complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel)>, 
		     std::vector< int >, double&);

  
  bool initializeHypothesisMap(const std::map<const std::string, bool>& hypMap);

  const double _pionMass;// = 0.13957;
  const double _etaMass;// = 0.547853;
  

};

#endif


