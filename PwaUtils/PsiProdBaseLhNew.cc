#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/PsiProdBaseLhNew.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


PsiProdBaseLhNew::PsiProdBaseLhNew(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, const std::vector<std::string>& hypVec, boost::shared_ptr<PsiToXGamStates> theStates) :
  AbsLhNew(theEvtList)
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false)
  ,_GammaEtaKey("GammaEta_")
  ,_GammaEta2Key("GammaEta2_")
  ,_GammaF0Key("GammaF0_")
  ,_GammaF1Key("GammaF1_")
  ,_GammaF2Key("GammaF2_")
  ,_EtaKey("Eta_")
  ,_Eta2Key("Eta2_")
  ,_F0Key("F0_")
  ,_F1Key("F1_")
  ,_F2Key("F2_")
  ,_phasespaceKey("Phasespace")
{
  initializeHypothesis();
}

PsiProdBaseLhNew::PsiProdBaseLhNew( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<PsiToXGamStates> theStates) :
  AbsLhNew(theLhPtr->getEventList())
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false)
  ,_GammaEtaKey("GammaEta_")
  ,_GammaEta2Key("GammaEta2_")
  ,_GammaF0Key("GammaF0_")
  ,_GammaF1Key("GammaF1_")
  ,_GammaF2Key("GammaF2_")
  ,_EtaKey("Eta_")
  ,_Eta2Key("Eta2_")
  ,_F0Key("F0_")
  ,_F1Key("F1_")
  ,_F2Key("F2_")
  ,_phasespaceKey("Phasespace")
{
  initializeHypothesis();  
}

PsiProdBaseLhNew::~PsiProdBaseLhNew()
{;
}

double PsiProdBaseLhNew::calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal){

  double result=0.;

  complex<double> JmpGmp(0.,0.);
  complex<double> JmpGmm(0.,0.);
  complex<double> JmmGmp(0.,0.);
  complex<double> JmmGmm(0.,0.);  

  //calculate all amplitudes
  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for (itMap=_hypMap.begin(); itMap!=_hypMap.end(); ++itMap){
    std::vector<std::string>::const_iterator itStr;

    for (itStr = itMap->second.begin(); itStr!= itMap->second.end(); ++itStr){
      boost::shared_ptr<AbsXdecAmp> currentDecAmp=_allAmpMap[*itStr];

      std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToXGamMag=theParamVal.MagLamLams[*itStr];
      std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToXGamPhi=theParamVal.PhiLamLams[*itStr];
      
      std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
      for (itMag=PsiToXGamMag.begin(); itMag!=PsiToXGamMag.end(); ++itMag){
	boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
        Spin helX=currentJPClamlam->lam1;

	complex<double> JpGpTmp(0.,0.);
	complex<double> JmGpTmp(0.,0.);
	complex<double> JpGmTmp(0.,0.);
	complex<double> JmGmTmp(0.,0.);
	complex<double> TmpDecAmpPos(0.,0.);
	complex<double> TmpDecAmpNeg(0.,0.);

	/*complex<double>*/ JpGpTmp = psiToXGammaAmp(1, helX, 1, theData, itMag->second, PsiToXGamPhi[currentJPClamlam]);
	/*complex<double>*/ JmGpTmp = psiToXGammaAmp(-1, helX, 1, theData, itMag->second, PsiToXGamPhi[currentJPClamlam]);
	
	/*complex<double>*/ JpGmTmp = currentJPClamlam->parityFactor*psiToXGammaAmp(1, -helX, -1, theData, itMag->second, PsiToXGamPhi[currentJPClamlam]);
	/*complex<double>*/ JmGmTmp = currentJPClamlam->parityFactor*psiToXGammaAmp(-1, -helX, -1, theData, itMag->second, PsiToXGamPhi[currentJPClamlam]);	

	/*complex<double>*/ TmpDecAmpPos = currentDecAmp->XdecAmp(helX, theData);
	/*complex<double>*/ TmpDecAmpNeg = TmpDecAmpPos;
        if (helX>0)  TmpDecAmpNeg = currentDecAmp->XdecAmp(-helX, theData);

	JmpGmp+=JpGpTmp*TmpDecAmpPos;
	JmpGmm+=JpGmTmp*TmpDecAmpNeg;
	JmmGmp+=JmGpTmp*TmpDecAmpPos;
	JmmGmm+=JmGmTmp*TmpDecAmpNeg;
      }
    }
  }

  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }
  return result;  

}

complex<double> PsiProdBaseLhNew::psiToXGammaAmp(Spin Minit, Spin lamX, Spin lamGamma, EvtDataNew* theData, 
                                                     double PsiToXGamMag, double PsiToXGamPhi ){
  complex<double> result(0.,0.);
   Spin lambda = lamX-lamGamma;
   complex<double> expiphiPsi(cos(PsiToXGamPhi), sin(PsiToXGamPhi)); 
   result = PsiToXGamMag*expiphiPsi*conj( theData->WignerDsProd[enumProdDfunc::Psi][1][Minit][lambda]  );

   return result;

}

void PsiProdBaseLhNew::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){ 

  if(_usePhasespace){
    fitVal.otherParams[_phasespaceKey]=0.02;
    fitErr.otherParams[_phasespaceKey]=0.015;
  }

  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for ( itMap = _hypMap.begin(); itMap != _hypMap.end(); ++itMap){
 
    std::vector< boost::shared_ptr<const JPClamlam> > currentlamLamAmps;   
    if (itMap->first==_GammaEtaKey) currentlamLamAmps=_theStatesPtr->PsiToEtaGammaLamLamStates();
    else if (itMap->first==_GammaEta2Key) currentlamLamAmps=_theStatesPtr->PsiToEta2GammaLamLamStates();
    else if (itMap->first==_GammaF0Key) currentlamLamAmps=_theStatesPtr->PsiToF0GammaLamLamStates();
    else if (itMap->first==_GammaF1Key) currentlamLamAmps=_theStatesPtr->PsiToF1GammaLamLamStates();
    else if (itMap->first==_GammaF2Key) currentlamLamAmps=_theStatesPtr->PsiToF2GammaLamLamStates();
    else {
      Alert << "Key\t" << itMap->first << "\t not supported!!!";
			   exit(1); 
    }

    getDefaultLamLamParams( itMap->second, currentlamLamAmps, fitVal, fitErr);   
    
 }

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::const_iterator it1;
  for (it1=_etaDecAmpMap.begin(); it1!=_etaDecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }    
  for (it1=_eta2DecAmpMap.begin(); it1!=_eta2DecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }
  for (it1=_f0DecAmpMap.begin(); it1!=_f0DecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }
  for (it1=_f1DecAmpMap.begin(); it1!=_f1DecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }
  for (it1=_f2DecAmpMap.begin(); it1!=_f2DecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }



}

void PsiProdBaseLhNew::print(std::ostream& os) const{

}


void  PsiProdBaseLhNew::initializeHypothesis(){

  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _GammaEtaKey.size(), _GammaEtaKey)== 0){
      Info << "eta hypothesis\t" << (*it) << "\t found" << endmsg;
      _GammaEtaHyps.push_back(*it);
      _hypMap[_GammaEtaKey].push_back(*it);
    }
    else if (it->compare(0, _GammaEta2Key.size(), _GammaEta2Key)== 0){
      Info << "eta2 hypothesis\t" << (*it) << "\t found" << endmsg;
      _GammaEta2Hyps.push_back(*it);
      _hypMap[_GammaEta2Key].push_back(*it);
    }
    else if (it->compare(0, _GammaF0Key.size(), _GammaF0Key)== 0){
      Info << "f0 hypothesis\t" << (*it) << "\t found" << endmsg; 
      _GammaF0Hyps.push_back(*it);
      _hypMap[_GammaF0Key].push_back(*it);
    }
    else if (it->compare(0, _GammaF1Key.size(), _GammaF1Key)== 0){
      Info << "f1 hypothesis\t" << (*it) << "\t found" << endmsg; 
      _GammaF1Hyps.push_back(*it);
      _hypMap[_GammaF1Key].push_back(*it);
    }
    else if (it->compare(0, _GammaF2Key.size(), _GammaF2Key)== 0){
      Info << "f2 hypothesis\t" << (*it) << "\t found" << endmsg; 
      _GammaF2Hyps.push_back(*it);
      _hypMap[_GammaF2Key].push_back(*it);
    }
    else if (it->compare(0, _phasespaceKey.size(), _phasespaceKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _usePhasespace=true;
    }
   
  } 

  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for ( itMap = _hypMap.begin(); itMap != _hypMap.end(); ++itMap){
    std::vector<std::string>::const_iterator itVec;
    for (itVec=itMap->second.begin(); itVec!=itMap->second.end(); ++itVec){
      std::cout << *itVec << std::endl;
    }
 }
}



void PsiProdBaseLhNew::getDefaultLamLamParams(const std::vector<std::string>& hyps, std::vector< boost::shared_ptr<const JPClamlam> > lamLamAmps, fitParamsNew& fitVal, fitParamsNew& fitErr){

  if (hyps.size()==0) return;

  std::vector<std::string>::const_iterator it;
  
  for (it=hyps.begin(); it!=hyps.end(); ++it){
    //create std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
    std::vector< boost::shared_ptr<const JPClamlam> >::const_iterator itLamLam;
    
    for(itLamLam=lamLamAmps.begin(); itLamLam!=lamLamAmps.end(); ++itLamLam){
      currentMagValMap[*itLamLam]=0.7;
      currentPhiValMap[*itLamLam]=0.;
      currentMagErrMap[*itLamLam]=0.6;
      currentPhiErrMap[*itLamLam]=0.3;
    }

    fitVal.MagLamLams[*it]=currentMagValMap;
    fitVal.PhiLamLams[*it]=currentPhiValMap;
    fitErr.MagLamLams[*it]=currentMagErrMap;
    fitErr.PhiLamLams[*it]=currentPhiErrMap;   
  }
  

}

void PsiProdBaseLhNew::checkParamVariation(fitParamsNew& theParamVal){
  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for (itMap=_hypMap.begin(); itMap!=_hypMap.end(); ++itMap){
    std::vector<std::string>::const_iterator itStr;

    for (itStr = itMap->second.begin(); itStr!= itMap->second.end(); ++itStr){
      boost::shared_ptr<AbsXdecAmp> currentDecAmp=_allAmpMap[*itStr];
      currentDecAmp->checkRecalculation(theParamVal);
    }  
  }  
  return;
}

void PsiProdBaseLhNew::cacheTheAmps(){
  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for (itMap=_hypMap.begin(); itMap!=_hypMap.end(); ++itMap){
    std::vector<std::string>::const_iterator itStr;

    for (itStr = itMap->second.begin(); itStr!= itMap->second.end(); ++itStr){
      boost::shared_ptr<AbsXdecAmp> currentDecAmp=_allAmpMap[*itStr];
      currentDecAmp->cacheAmplitudes();
    }  
  }  
  return;
}

void PsiProdBaseLhNew::updateFitParams(fitParamsNew& theParamVal){

  std::map< std::string,std::vector<std::string> >::const_iterator itMap;
  for (itMap=_hypMap.begin(); itMap!=_hypMap.end(); ++itMap){
    std::vector<std::string>::const_iterator itStr;
    
    for (itStr = itMap->second.begin(); itStr!= itMap->second.end(); ++itStr){
      boost::shared_ptr<AbsXdecAmp> currentDecAmp=_allAmpMap[*itStr];
      currentDecAmp->updateFitParams(theParamVal);
    }
  }
}
