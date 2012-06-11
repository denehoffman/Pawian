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
  ,_GammaF1Key("GammaF1_")
  ,_EtaKey("Eta_")
  ,_Eta2Key("Eta2_")
  ,_F1Key("F1_")
  ,_phasespaceKey("Phasespace")
{
  initializeHypothesis();
 
}

PsiProdBaseLhNew::PsiProdBaseLhNew( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<PsiToXGamStates> theStates ) :
  AbsLhNew(theLhPtr->getEventList())
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false)
  ,_GammaEtaKey("GammaEta_")
  ,_GammaEta2Key("GammaEta2_")
  ,_GammaF1Key("GammaF1_")
  ,_EtaKey("Eta_")
  ,_Eta2Key("Eta2_")
  ,_F1Key("F1_")
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

  //calculate gamma eta amplitudes

  std::vector<std::string>::const_iterator itStr;
  for (itStr= _GammaEtaHyps.begin(); itStr!= _GammaEtaHyps.end(); ++itStr){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEtaGamMag=theParamVal.MagLamLams[*itStr];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEtaGamPhi=theParamVal.PhiLamLams[*itStr];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=PsiToEtaGamMag.begin(); itMag!=PsiToEtaGamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiToEtaGamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiToEtaGamPhi[currentJPClamlam];      
    }

    complex<double>  JpGpTmp = ParityProdMap[0][1]*psiToXGammaAmp(1, 0, 0, 1, theData, MagProdMap[0][1], PhiProdMap[0][1]);
    complex<double>  JpGmTmp = ParityProdMap[0][-1]*psiToXGammaAmp(1, 0, 0, -1, theData, MagProdMap[0][-1], PhiProdMap[0][-1]);
    complex<double>  JmGpTmp = ParityProdMap[0][1]*psiToXGammaAmp(-1, 0, 0, 1, theData, MagProdMap[0][1], PhiProdMap[0][1]);
    complex<double>  JmGmTmp = ParityProdMap[0][-1]*psiToXGammaAmp(-1, 0, 0, -1, theData, MagProdMap[0][-1], PhiProdMap[0][-1]);

    
    boost::shared_ptr<AbsXdecAmp> currentEtaDecAmp=_etaDecAmpMap[*itStr];
    complex<double>  TmpDecAmp=currentEtaDecAmp->XdecAmp(0, theData, theParamVal);    
 
    JmpGmp+=JpGpTmp*TmpDecAmp;
    JmpGmm+=JpGmTmp*TmpDecAmp;
    JmmGmp+=JmGpTmp*TmpDecAmp;
    JmmGmm+=JmGmTmp*TmpDecAmp;
  }
  
  //calculate gamma f1 amplitudes
  for (itStr= _GammaF1Hyps.begin(); itStr!= _GammaF1Hyps.end(); ++itStr){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToF1GamMag=theParamVal.MagLamLams[*itStr];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToF1GamPhi=theParamVal.PhiLamLams[*itStr];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=PsiToF1GamMag.begin(); itMag!=PsiToF1GamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiToF1GamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiToF1GamPhi[currentJPClamlam];      
    }

    std::map<Spin,complex<double> > JpGpTmpMap;
    std::map<Spin,complex<double> > JpGmTmpMap;
    std::map<Spin,complex<double> > JmGpTmpMap;
    std::map<Spin,complex<double> > JmGmTmpMap;
    std::map<Spin,complex<double> > TmpDecAmp;

    boost::shared_ptr<AbsXdecAmp> currentF1DecAmp=_f1DecAmpMap[*itStr];
    
    for (Spin helf1=0; helf1<2; helf1++){
      JpGpTmpMap[helf1]= ParityProdMap[helf1][1]*psiToXGammaAmp(1, 1, helf1, 1, theData, MagProdMap[helf1][1], PhiProdMap[helf1][1]);
      JmGpTmpMap[helf1]= ParityProdMap[helf1][1]*psiToXGammaAmp(-1, 1, helf1, 1, theData, MagProdMap[helf1][1], PhiProdMap[helf1][1]);
      TmpDecAmp[helf1] = currentF1DecAmp->XdecAmp(helf1, theData, theParamVal);

      if(helf1>0){
        JpGmTmpMap[-helf1]= ParityProdMap[-helf1][-1]*psiToXGammaAmp(1, 1, -helf1, -1, theData, MagProdMap[-helf1][-1], PhiProdMap[-helf1][-1]);
        JmGmTmpMap[-helf1]= ParityProdMap[-helf1][-1]*psiToXGammaAmp(-1, 1, -helf1, -1, theData, MagProdMap[-helf1][-1], PhiProdMap[-helf1][-1]);
        TmpDecAmp[-helf1] = currentF1DecAmp->XdecAmp(-helf1, theData, theParamVal);
      }
    }

    for (Spin helf1=-1; helf1<2; helf1++){
      JmpGmp+=JpGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmpGmm+=JpGmTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmp+=JmGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmm+=JmGmTmpMap[helf1]*TmpDecAmp[helf1];

    }

  }



  //calculate gamma eta2 amplitudes
  for (itStr= _GammaEta2Hyps.begin(); itStr!= _GammaEta2Hyps.end(); ++itStr){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamMag=theParamVal.MagLamLams[*itStr];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamPhi=theParamVal.PhiLamLams[*itStr];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=PsiToEta2GamMag.begin(); itMag!=PsiToEta2GamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiToEta2GamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiToEta2GamPhi[currentJPClamlam];      
    }

    std::map<Spin,complex<double> > JpGpTmpMap;
    std::map<Spin,complex<double> > JpGmTmpMap;
    std::map<Spin,complex<double> > JmGpTmpMap;
    std::map<Spin,complex<double> > JmGmTmpMap;
    std::map<Spin,complex<double> > TmpDecAmp;

    boost::shared_ptr<AbsXdecAmp> currentEta2DecAmp=_eta2DecAmpMap[*itStr];
    
    for (Spin heleta2=0; heleta2<3; heleta2++){
      JpGpTmpMap[heleta2]= ParityProdMap[heleta2][1]*psiToXGammaAmp(1, 2, heleta2, 1, theData, MagProdMap[heleta2][1], PhiProdMap[heleta2][1]);
      JmGpTmpMap[heleta2]= ParityProdMap[heleta2][1]*psiToXGammaAmp(-1, 2, heleta2, 1, theData, MagProdMap[heleta2][1], PhiProdMap[heleta2][1]);
      TmpDecAmp[heleta2] = currentEta2DecAmp->XdecAmp(heleta2, theData, theParamVal);

      if(heleta2>0){
        JpGmTmpMap[-heleta2]= ParityProdMap[-heleta2][-1]*psiToXGammaAmp(1, 2, -heleta2, -1, theData, MagProdMap[-heleta2][-1], PhiProdMap[-heleta2][-1]);
        JmGmTmpMap[-heleta2]= ParityProdMap[-heleta2][-1]*psiToXGammaAmp(-1, 2, -heleta2, -1, theData, MagProdMap[-heleta2][-1], PhiProdMap[-heleta2][-1]);
        TmpDecAmp[-heleta2] = currentEta2DecAmp->XdecAmp(-heleta2, theData, theParamVal);
      }
    }

    for (Spin heleta2=-2; heleta2<3; heleta2++){
      JmpGmp+=JpGpTmpMap[heleta2]*TmpDecAmp[heleta2];
      JmpGmm+=JpGmTmpMap[heleta2]*TmpDecAmp[heleta2];
      JmmGmp+=JmGpTmpMap[heleta2]*TmpDecAmp[heleta2];
      JmmGmm+=JmGmTmpMap[heleta2]*TmpDecAmp[heleta2];

    }

  }



  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);


  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }

//   Info << "result:\t" << result << endmsg;
  return result;  

}

complex<double> PsiProdBaseLhNew::psiToXGammaAmp(Spin Minit, Spin jX, Spin lamX, Spin lamGamma, EvtDataNew* theData, 
                                                     double PsiToXGamMag, double PsiToXGamPhi ){
   Spin lambda = lamX-lamGamma;
   complex<double> expiphiPsi(cos(PsiToXGamPhi), sin(PsiToXGamPhi));
   complex<double> result = PsiToXGamMag*expiphiPsi*conj( theData->WignerDsProd[enumProdDfunc::Psi][1][Minit][lambda]  );

   return result;

}

void PsiProdBaseLhNew::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){ 

  if(_usePhasespace){
    fitVal.otherParams[_phasespaceKey]=0.02;
    fitErr.otherParams[_phasespaceKey]=0.015;
  }

  getDefaultLamLamParams(_GammaEtaHyps, _theStatesPtr->PsiToEtaGammaLamLamStates(), fitVal, fitErr);  
  getDefaultLamLamParams(_GammaEta2Hyps, _theStatesPtr->PsiToEta2GammaLamLamStates(), fitVal, fitErr); 
  getDefaultLamLamParams(_GammaF1Hyps, _theStatesPtr->PsiToF1GammaLamLamStates(), fitVal, fitErr);

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::const_iterator it1;
  for (it1=_etaDecAmpMap.begin(); it1!=_etaDecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }    

  for (it1=_f1DecAmpMap.begin(); it1!=_f1DecAmpMap.end(); ++it1){
    it1->second->getDefaultParams(fitVal, fitErr);
  }

  for (it1=_eta2DecAmpMap.begin(); it1!=_eta2DecAmpMap.end(); ++it1){
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

    }
    else if (it->compare(0, _GammaF1Key.size(), _GammaF1Key)== 0){
      Info << "f1 hypothesis\t" << (*it) << "\t found" << endmsg;
      _GammaF1Hyps.push_back(*it);
    }
    else if (it->compare(0, _GammaEta2Key.size(), _GammaEta2Key)== 0){
      Info << "eta2 hypothesis\t" << (*it) << "\t found" << endmsg;
      _GammaEta2Hyps.push_back(*it);
    }
    else if (it->compare(0, _phasespaceKey.size(), _phasespaceKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _usePhasespace=true;
    }
   
  } 

}



void PsiProdBaseLhNew::getDefaultLamLamParams(std::vector<std::string>& hyps, std::vector< boost::shared_ptr<const JPClamlam> > lamLamAmps, fitParamsNew& fitVal, fitParamsNew& fitErr){

  if (hyps.size()==0) return;

  std::vector<std::string>::iterator it;
  
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
