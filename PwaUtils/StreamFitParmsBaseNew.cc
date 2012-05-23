#include <fstream>

#include "PwaUtils/StreamFitParmsBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"
#include "ErrLogger/ErrLogger.hh"

StreamFitParmsBaseNew::StreamFitParmsBaseNew(std::string& filePath, boost::shared_ptr<AbsLhNew> theLhPtr) :
  AbsFitParamStreamer(filePath)
{
  theLhPtr->getDefaultParams(_paramVal,_paramErr);
  fillParams();

}

StreamFitParmsBaseNew::~StreamFitParmsBaseNew(){;}

void StreamFitParmsBaseNew::fillParams(){
  const std::string magSuffix="Mag";
  const std::string phiSuffix="Phi"; 
  const std::string massSuffix="Mass";
  const std::string widthSuffix="Width";
  const std::string gFactorSuffix="gFactor";

  fillLamLamAmps(_paramVal.MagLamLams, _paramErr.MagLamLams, magSuffix);
  fillLamLamAmps(_paramVal.PhiLamLams, _paramErr.PhiLamLams, phiSuffix);
  fillLSAmps(_paramVal.Mags, _paramErr.Mags, magSuffix);
  fillLSAmps(_paramVal.Phis, _paramErr.Phis, phiSuffix);
  fillDoubles(_paramVal.Masses, _paramErr.Masses, massSuffix);
  fillDoubles(_paramVal.Widths, _paramErr.Widths, widthSuffix);
  fillDoubles(_paramVal.gFactors, _paramErr.gFactors, gFactorSuffix);
}

void StreamFitParmsBaseNew::fillLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap, const std::string& suffix){

  mapStrJPCLamLam::iterator itLamLamMap;
  for( itLamLamMap=valMap.begin(); itLamLamMap!=valMap.end(); ++itLamLamMap){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLam;

    for ( itLamLam=itLamLamMap->second.begin(); itLamLam!=itLamLamMap->second.end();  ++itLamLam){
     

      std::string theKey=itLamLam->first->name()+itLamLamMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;
      
      stringPairIter=_stringPairMap.find(theKey);
      
      if ( stringPairIter != _stringPairMap.end() ){
	Info << "key\t" << theKey << "\tfound" << endmsg;
	double val=stringPairIter->second.first;
	Info << "replace val by " << val << endmsg;
	double err=stringPairIter->second.second;
	Info << "replace err by " << err << endmsg;

	valMap[itLamLamMap->first][itLamLam->first] = val;
	errMap[itLamLamMap->first][itLamLam->first] = err;

      }
    }
  }
}

void StreamFitParmsBaseNew::fillLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap, const std::string& suffix){

  mapStrJPCLS::iterator itLSMap;
  for( itLSMap=valMap.begin(); itLSMap!=valMap.end(); ++itLSMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itLS;

    for ( itLS=itLSMap->second.begin(); itLS!=itLSMap->second.end();  ++itLS){
     

      std::string theKey=itLS->first->name()+itLSMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;
      
      stringPairIter=_stringPairMap.find(theKey);
      
      if ( stringPairIter != _stringPairMap.end() ){
	Info << "key\t" << theKey << "\tfound" << endmsg;
	double val=stringPairIter->second.first;
	Info << "replace val by " << val << endmsg;
	double err=stringPairIter->second.second;
	Info << "replace err by " << err << endmsg;

	valMap[itLSMap->first][itLS->first] = val;
	errMap[itLSMap->first][itLS->first] = err;

      }
    }
  }
}

void StreamFitParmsBaseNew::fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap, const std::string& suffix){

  mapStrDouble::iterator itMap;
  for( itMap=valMap.begin(); itMap!=valMap.end(); ++itMap){

    std::string theKey=itMap->first+suffix;
    Info << "theKey=\t" << theKey << endmsg;
    StringPairMap::const_iterator stringPairIter;
    
    stringPairIter=_stringPairMap.find(theKey);

    if ( stringPairIter != _stringPairMap.end() ){
      Info << "key\t" << theKey << "\tfound" << endmsg;
      double val=stringPairIter->second.first;
      Info << "replace val by " << val << endmsg;
      double err=stringPairIter->second.second;
      Info << "replace err by " << err << endmsg;
      
      valMap[itMap->first] = val;
      errMap[itMap->first] = err;
      
    }
    
  }

}

void StreamFitParmsBaseNew::fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index){

  StringPairMap::const_iterator stringPairIter=_stringPairMap.find(suffix);
  if ( stringPairIter != _stringPairMap.end() ){
    theValMap[index]=stringPairIter->second.first;
    theErrMap[index]=stringPairIter->second.second;    
  }

}
