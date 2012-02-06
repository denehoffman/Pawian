#ifndef _JpsiGamKsKlKKProdLh_H
#define _JpsiGamKsKlKKProdLh_H

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
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKDynamicModelParams.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"




class JpsiGamKsKlKKProdLh : public AbsLh{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamKsKlKKProdLh(boost::shared_ptr<const EvtDataBaseList>, const std::map<const std::string, bool>& hypMap);
  JpsiGamKsKlKKProdLh(boost::shared_ptr<AbsLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~JpsiGamKsKlKKProdLh();

  virtual AbsLh* clone_() const {
    return new JpsiGamKsKlKKProdLh(_evtListPtr, _hypMap);
  }
  

  
  virtual double calcEvtIntensity( EvtData* theData, fitParams& theParamVal);
  
  //Getters:
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;
  
  void useCommonProductionPhase( bool commonPhase ){_useCommonProductionPhase=commonPhase;}
  void massIndependentFit( bool massIndep){ _massIndependentFit=massIndep;  }
  
protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, fitParams&  theParamVal, EvtData* theData);  

  
  virtual complex<double> etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, EvtData* theData, 
				      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaMag, 
				      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaPhi );
  
  virtual complex<double> f0GammaAmp(Spin Minit, Spin Mgamma, EvtData* theData, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdMag, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdPhi,
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag,  
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi,
				     bool useCommonProductionPhase, dynamicModelParams& dynModPars );
  
  virtual complex<double> etaToPhiPhiTo4KAmp(EvtData* theData);
  
  virtual complex<double> f0ToPhiPhiTo4KAmp( EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag, 
					     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi, dynamicModelParams& dynModPars    );
  
  
  virtual complex<double> f2GammaAmp(Spin Minit, Spin Mgamma, Spin fJSpin, EvtData* theData, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdMag, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdPhi,
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecMag,  
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecPhi,
				     double mass, double width, bool useCommonProductionPhase, dynamicModelParams& dynModPars );
  complex<double> f2ToPhiPhiTo4KAmp( EvtData* theData, Spin f2Lambda, 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecMag , 
				     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecPhi,
				     dynamicModelParams& dynModPars);
  

  complex<double> phiphiTo4KAmp(  EvtData* theData, Spin lambdaPhi1, Spin lambdaPhi2  );

  bool  initializeFlatteModel( dynamicModelParams &theDynModel, const Vector4<double> &fv2Phi, double mass, double gPhiPhi, double gKK );
  bool  initializeBreitWignerModel(dynamicModelParams &theDynModel,const Vector4<double> &fv2Phi, double mass, double width, dynamicModelParams::enumDynamicModel theModel  );
  
  bool _etacHyp;
  bool _eta2225Hyp;
  bool _f02020Hyp;
  bool _f02020FlatteHyp;
  
  bool _f22300Hyp;
  bool _eta21870Hyp;
  bool _f1Hyp;
  
  bool _usePhasespace;
  bool _useCommonProductionPhase;
  bool _massIndependentFit;
  
  std::map<const std::string, bool> _hypMap;

 
  

private:
  
  bool initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   );

  const double _phiMass;// = 1.019455;
  const double _kaonMass;// = 0.493677;
  

};

#endif
