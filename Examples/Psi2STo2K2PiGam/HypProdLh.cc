#include <getopt.h>
#include <fstream>

#include "Examples/Psi2STo2K2PiGam/HypProdLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

HypProdLh::HypProdLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList) :
  AbsPsi2STo2K2PiGamLh(theEvtList)
  ,_nFitParams(0) 
{
  setUp(); 
}

HypProdLh::HypProdLh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr) :
  AbsPsi2STo2K2PiGamLh(theLhPtr->getEventList())
  ,_nFitParams(0) 
{
  setUp();
}

HypProdLh::~HypProdLh()
{;
}

bool  HypProdLh::equalChic0DecParams(){
  return true;
}



complex<double> HypProdLh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result(1.,0.);

  return result; 

}



void HypProdLh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  _fitParams2K2PiGam.setMnUsrParamsFlatteMass(upar, startVal, errVal, paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace));

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }



}



int HypProdLh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  int counter=0;

  counter=_fitParams2K2PiGam.setFitParamFlatteMass(theParamVal, par, counter, paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace));

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, (*itAmps));
  }

  return counter;
}

unsigned int  HypProdLh::nFitParams(){

  return _nFitParams;
}


void HypProdLh::print(std::ostream& os) const{
  os << "HypProdLh::print\n";
}

void HypProdLh::printCurrentFitResult(param2K2PiGam& theParamVal){

  DebugMsg<< "phaseSpace: " << theParamVal.phaseSpace << endmsg;


  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< paramEnum2K2PiGam::name(*itAmps) << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }  
  }

}

void HypProdLh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  std::string tmpStringPs="phaseSpace"+suffix;
  os << tmpStringPs << "\t" << theParamVal.phaseSpace  << "\t" << 0. << std::endl;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];

      std::string tmpStringDec=(*itJPCLS)->name()+paramEnum2K2PiGam::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }

}


void HypProdLh::setUp(){

  _ampVec.push_back(paramEnum2K2PiGam::ChiGam);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  _nFitParams+=1; //phase space

}
