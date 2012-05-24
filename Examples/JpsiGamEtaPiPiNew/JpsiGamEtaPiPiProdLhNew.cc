#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiProdLhNew.hh"
#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


JpsiGamEtaPiPiProdLhNew::JpsiGamEtaPiPiProdLhNew(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates) :
  PsiProdBaseLhNew(theEvtList, hypVec, theStates)
  ,_jpsiGamEtaPiPiStatesPtr(theStates)
{
  initializeHypothesis();
 
}

JpsiGamEtaPiPiProdLhNew::JpsiGamEtaPiPiProdLhNew( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates ) :
  PsiProdBaseLhNew(theLhPtr->getEventList(), hypVec, theStates)
  ,_jpsiGamEtaPiPiStatesPtr(theStates)
{  
  initializeHypothesis();  
}

JpsiGamEtaPiPiProdLhNew::~JpsiGamEtaPiPiProdLhNew()
{;
}

void JpsiGamEtaPiPiProdLhNew::print(std::ostream& os) const{

}

void  JpsiGamEtaPiPiProdLhNew::initializeHypothesis(){

  std::vector<std::string>::const_iterator it;

  for (it=_GammaEtaHyps.begin(); it!=_GammaEtaHyps.end(); ++it){
    size_t pos=it->find(_EtaKey);
    std::string etaDecAmpName=it->substr(pos);
    _etaDecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XDecAmpBase( etaDecAmpName, _hypVec, _jpsiGamEtaPiPiStatesPtr, Spin(0)) );
  }

  for (it=_GammaF1Hyps.begin(); it!=_GammaF1Hyps.end(); ++it){
    size_t pos=it->find(_F1Key);
    std::string f1DecAmpName=it->substr(pos);
    _f1DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XDecAmpBase( f1DecAmpName, _hypVec, _jpsiGamEtaPiPiStatesPtr, Spin(1)) );
  }

  for (it=_GammaEta2Hyps.begin(); it!=_GammaEta2Hyps.end(); ++it){
    size_t pos=it->find(_Eta2Key);
    std::string eta2DecAmpName=it->substr(pos);
    _eta2DecAmpMap[*it]=boost::shared_ptr<AbsXdecAmp>(new XDecAmpBase( eta2DecAmpName, _hypVec, _jpsiGamEtaPiPiStatesPtr, Spin(2)) );
  }

//   std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::const_iterator it1;
//   for (it1=_etaDecAmpMap.begin(); it1!=_etaDecAmpMap.end(); ++it1){
//     Info << "eta dec amp:\t" << it1->second->name() << endmsg;
//   }
//   for (it1=_f1DecAmpMap.begin(); it1!=_f1DecAmpMap.end(); ++it1){
//     Info << "f1 dec amp:\t" << it1->second->name() << endmsg;
//   }
//   for (it1=_eta2DecAmpMap.begin(); it1!=_eta2DecAmpMap.end(); ++it1){
//     Info << "eta2 dec amp:\t" << it1->second->name() << endmsg;
//   }
}

