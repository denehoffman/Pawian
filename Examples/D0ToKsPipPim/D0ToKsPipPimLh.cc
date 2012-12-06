#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToKsPipPimLh.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"
#include "Examples/D0ToKsPipPim/D0ToPiPiSKDec.hh"
#include "Examples/D0ToKsPipPim/D0ToKPiSPiDec.hh"

#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


D0ToKsPipPimLh::D0ToKsPipPimLh(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsLhNew(theEvtList)
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false)
  ,_phasespaceKey("Phasespace") 
  ,_pipiSHyp(false)
  ,_pipiSKey("PiPiSWave")
  ,_KpiSHyp(false)
  ,_KpiSKey("KPiSWave")
{
  initializeHypothesis();
}

D0ToKsPipPimLh::D0ToKsPipPimLh( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsLhNew(theLhPtr->getEventList())
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false)
  ,_phasespaceKey("Phasespace") 
  ,_pipiSHyp(false)
  ,_pipiSKey("PiPiSWave")
  ,_KpiSHyp(false)
  ,_KpiSKey("KPiSWave")
{
  initializeHypothesis();  
}

D0ToKsPipPimLh::~D0ToKsPipPimLh()
{;
}

double D0ToKsPipPimLh::calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal){

  double result=0.;

  complex<double> tmpAmp(0.,0.);

 //calculate all amplitudes
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it = _allDecAmpMap.begin(); it != _allDecAmpMap.end(); ++it){
    boost::shared_ptr<AbsXdecAmp> currentDecAmp=it->second;
  
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMag=theParamVal.Mags[it->first];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > thePhi=theParamVal.Phis[it->first];
 
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag; 
 
    for ( itXMag=theMag.begin(); itXMag!=theMag.end(); ++itXMag ){
      boost::shared_ptr<const JPCLS> XState=itXMag->first;
      double theXMag=itXMag->second;
      double theXPhi=thePhi[XState];

// #ifdef _OPENMP
// #pragma omp critical
//   {
// #endif 
//       std::cout << "theXMag:\t" << theXMag <<"\ttheXPhi:\t" << theXPhi << std::endl;
// #ifdef _OPENMP
//   }
// #endif
     
      complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
      
      complex<double> amp = theXMag*expiphiX*currentDecAmp->XdecAmp(0, theData);
      tmpAmp+= amp;
    }
  }
  
  result=norm(tmpAmp);
  
  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];  
  }
  
  return result;  

}

void D0ToKsPipPimLh::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){ 

  if(_usePhasespace){
    fitVal.otherParams[_phasespaceKey]=0.02;
    fitErr.otherParams[_phasespaceKey]=0.015;
  }

  std::map< std::string, std::vector<std::string> >::const_iterator itMap;
  for ( itMap = _hypMaps.begin(); itMap != _hypMaps.end(); ++itMap){

  
    std::vector< boost::shared_ptr<const JPCLS> > currentLSAmps;
 
    if (itMap->first == _pipiSKey){
      currentLSAmps=_theStatesPtr->D0Tof0KStates();
    }
    else if (itMap->first == _KpiSKey){
      currentLSAmps=_theStatesPtr->D0ToK0KStates();
    }
    else {
      Alert << "Key\t" << itMap->first << "\t not supported!!!";
                           exit(1); 
    }

    getDefaultLSParams(itMap->second, currentLSAmps, fitVal, fitErr);

    std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::const_iterator itDecAmp;
    for (itDecAmp=_allDecAmpMap.begin(); itDecAmp!=_allDecAmpMap.end(); ++itDecAmp){
      itDecAmp->second->getDefaultParams(fitVal, fitErr);
    }    
  }

}

void D0ToKsPipPimLh::print(std::ostream& os) const{

}

void D0ToKsPipPimLh::getDefaultLSParams(const std::vector<std::string>& hyps, 
					std::vector< boost::shared_ptr<const JPCLS> > lsAmps, 
					fitParamsNew& fitVal, fitParamsNew& fitErr){  

  if (hyps.size()==0) return;

  std::vector<std::string>::const_iterator it;
  
  for (it=hyps.begin(); it!=hyps.end(); ++it){


    //create std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    for(itLS=lsAmps.begin(); itLS!=lsAmps.end(); ++itLS){
      currentMagValMap[*itLS]=0.7;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.6;
      currentPhiErrMap[*itLS]=0.3;
    }

    fitVal.Mags[*it]=currentMagValMap;
    fitVal.Phis[*it]=currentPhiValMap;
    fitErr.Mags[*it]=currentMagErrMap;
    fitErr.Phis[*it]=currentPhiErrMap;   
  }

}

void  D0ToKsPipPimLh::initializeHypothesis(){
  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _pipiSKey.size(), _pipiSKey)== 0){
      if(_pipiSHyp){
	Alert << "pipiS wave already initialized. 2 pipiS wave are forbidden!!!" << endmsg;
	exit(0);
      }
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _hypMaps[_pipiSKey].push_back(*it);
      _allDecAmpMap[(*it)]=boost::shared_ptr<AbsXdecAmp>(new D0ToPiPiSKDec( (*it), _hypVec, _theStatesPtr));
      _pipiSHyp=true;
    }
    else if (it->compare(0, _KpiSKey.size(), _KpiSKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _hypMaps[_KpiSKey].push_back(*it);
      _allDecAmpMap[(*it)]=boost::shared_ptr<AbsXdecAmp>(new D0ToKPiSPiDec( (*it), _hypVec, _theStatesPtr));
      _KpiSHyp=true;
    }
    else if (it->compare(0, _phasespaceKey.size(), _phasespaceKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _usePhasespace=true;
    }

  }

}




