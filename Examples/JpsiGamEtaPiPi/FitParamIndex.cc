#include "Examples/JpsiGamEtaPiPi/FitParamIndex.hh"


FitParamIndex::FitParamIndex(fitParams& theParams){
  
  int indexCounter=0;
  
  //iterate amplitude magnitueds
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& startMagMap=theParams.Mags;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::map jpclsMap = itMagMap->second;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator jpclsMapIter;
    for (jpclsMapIter=jpclsMap.begin(); jpclsMapIter!=jpclsMap.end(); ++jpclsMapIter){
      _Mags[itMagMap->first][jpclsMapIter->first] = indexCounter;
      indexCounter++;
    }
  }
  
  //loop again and fill phi index
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::map jpclsMap = itMagMap->second;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator jpclsMapIter;
    for (jpclsMapIter=jpclsMap.begin(); jpclsMapIter!=jpclsMap.end(); ++jpclsMapIter){
      _Phis[itMagMap->first][jpclsMapIter->first] = indexCounter;
      indexCounter++;
    }
  }
  
}
