#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsLhNew.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

AbsLhNew::AbsLhNew(boost::shared_ptr<const EvtDataBaseListNew> theEvtList) :
  _evtListPtr(theEvtList)
  ,_cacheAmps(false)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLhNew::AbsLhNew(boost::shared_ptr<AbsLhNew> theAbsLhPtr):
  _evtListPtr(theAbsLhPtr->getEventList())
  ,_cacheAmps(false)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLhNew::~AbsLhNew()
{
}

double AbsLhNew::calcLogLh(fitParamsNew& theParamVal){

  if (_cacheAmps) checkParamVariation(theParamVal); 
  updateFitParams(theParamVal);
  
  double logLH=0.;
  double logLH_data=0.;
  double weightSum=0.;
  double LH_mc=0.;


#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i=0; i<_evtDataVec.size(); ++i){
    EvtDataNew* currentEvtData=_evtDataVec[i];
    double intensity=calcEvtIntensity(currentEvtData, theParamVal);

#ifdef _OPENMP
#pragma omp critical
      {
#endif
    if (intensity>0.) logLH_data+=(currentEvtData->evtWeight)*log(intensity);
    weightSum+= currentEvtData->evtWeight;

#ifdef _OPENMP
       }
#endif
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i=0; i<_evtMCVec.size(); ++i){
    EvtDataNew* currentEvtData=_evtMCVec[i];
    double intensity=calcEvtIntensity(currentEvtData, theParamVal);

#ifdef _OPENMP
#pragma omp critical
      {
#endif
        LH_mc+=intensity;
#ifdef _OPENMP
       }
#endif
  }


  double logLH_mc_Norm=0.;  
  if (LH_mc>0.) logLH_mc_Norm=log(LH_mc/_evtMCVec.size());
  logLH=0.5*weightSum *(LH_mc/_evtMCVec.size()-1.)*(LH_mc/_evtMCVec.size()-1.)
    -logLH_data
    +weightSum*logLH_mc_Norm;
  
  Info << "current LH = " << logLH << endmsg;
  return logLH;
  
}

void AbsLhNew::setHyps( const std::map<const std::string, bool>& theMap, bool& theHyp, std::string& theKey){

  std::map<const std::string, bool>::const_iterator iter= theMap.find(theKey);
  
  if (iter !=theMap.end()){
    theHyp= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << theHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << theKey << " does not exist!!!" <<endmsg;
    exit(0);
  }
}


void AbsLhNew::cacheAmplitudes(){
  _cacheAmps=true;
  cacheTheAmps();
}


void AbsLhNew::checkParamVariation(fitParamsNew& theParamVal){

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it = _allDecAmpMap.begin(); it != _allDecAmpMap.end(); ++it){
    it->second->checkRecalculation(theParamVal);
  }

  return;
}

void AbsLhNew::cacheTheAmps(){

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it = _allDecAmpMap.begin(); it != _allDecAmpMap.end(); ++it){
    it->second->cacheAmplitudes();
  }
  return;
}

void AbsLhNew::updateFitParams(fitParamsNew& theParamVal){
  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it = _allDecAmpMap.begin(); it != _allDecAmpMap.end(); ++it){
    it->second->updateFitParams(theParamVal);
  }
  return;
}
