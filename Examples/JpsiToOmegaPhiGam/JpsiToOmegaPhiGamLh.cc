#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiToOmegaPhiGam/JpsiToOmegaPhiGamLh.hh"
#include "Examples/JpsiToOmegaPhiGam/XToOmegaPhiDecAmps.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


JpsiToOmegaPhiGamLh::JpsiToOmegaPhiGamLh(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToOmegaPhiGamStates> theStates) :
  PsiProdBaseLhNew(theEvtList, hypVec, theStates)
  ,_jpsiToOmegaPhiGamStatesPtr(theStates)
{
  initializeHypothesis();
 
}

JpsiToOmegaPhiGamLh::JpsiToOmegaPhiGamLh( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToOmegaPhiGamStates> theStates) :
  PsiProdBaseLhNew(theLhPtr->getEventList(), hypVec, theStates)
  ,_jpsiToOmegaPhiGamStatesPtr(theStates)
{  
  initializeHypothesis();  
}

JpsiToOmegaPhiGamLh::~JpsiToOmegaPhiGamLh()
{;
}

void JpsiToOmegaPhiGamLh::print(std::ostream& os) const{

}

void  JpsiToOmegaPhiGamLh::initializeHypothesis(){

  std::vector<std::string>::const_iterator it;

  for (it=_GammaEtaHyps.begin(); it!=_GammaEtaHyps.end(); ++it){
    size_t pos=it->find(_EtaKey);
    std::string etaDecAmpName=it->substr(pos);
    _etaDecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XToOmegaPhiDecAmps( etaDecAmpName, _hypVec, _jpsiToOmegaPhiGamStatesPtr, Spin(0), -1 ));
    _allAmpMap[*it]=_etaDecAmpMap[*it];
  }
  for (it=_GammaEta2Hyps.begin(); it!=_GammaEta2Hyps.end(); ++it){
    size_t pos=it->find(_Eta2Key);
    std::string eta2DecAmpName=it->substr(pos);
    _eta2DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XToOmegaPhiDecAmps( eta2DecAmpName, _hypVec, _jpsiToOmegaPhiGamStatesPtr, Spin(2), -1 ));
    _allAmpMap[*it]=_eta2DecAmpMap[*it];
  }

  for (it=_GammaF0Hyps.begin(); it!=_GammaF0Hyps.end(); ++it){
    size_t pos=it->find(_F0Key);
    std::string f0DecAmpName=it->substr(pos);
    _f0DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XToOmegaPhiDecAmps( f0DecAmpName, _hypVec, _jpsiToOmegaPhiGamStatesPtr, Spin(0), 1 ));
    _allAmpMap[*it]=_f0DecAmpMap[*it];
  }

  for (it=_GammaF1Hyps.begin(); it!=_GammaF1Hyps.end(); ++it){
    size_t pos=it->find(_F1Key);
    std::string f1DecAmpName=it->substr(pos);
    _f1DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XToOmegaPhiDecAmps( f1DecAmpName, _hypVec, _jpsiToOmegaPhiGamStatesPtr, Spin(1), 1 ));
     _allAmpMap[*it]=_f1DecAmpMap[*it];
  }

  for (it=_GammaF2Hyps.begin(); it!=_GammaF2Hyps.end(); ++it){
    size_t pos=it->find(_F2Key);
    std::string f2DecAmpName=it->substr(pos);
    _f2DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XToOmegaPhiDecAmps( f2DecAmpName, _hypVec, _jpsiToOmegaPhiGamStatesPtr, Spin(2),1  ));
    _allAmpMap[*it]=_f2DecAmpMap[*it];
  }

}



