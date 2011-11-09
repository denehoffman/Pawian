#include <fstream>

#include "Examples/Psi2SToKpKmPiGam/StreamChic1ToKpKmPiGamFitParms.hh"
#include "ErrLogger/ErrLogger.hh"

StreamChic1ToKpKmPiGamFitParms::StreamChic1ToKpKmPiGamFitParms(std::string& filePath) :
  AbsFitParamStreamer(filePath),
  _fitParamsKpKmPiGam()
{
  fillParams();

}

StreamChic1ToKpKmPiGamFitParms::~StreamChic1ToKpKmPiGamFitParms(){;}

void StreamChic1ToKpKmPiGamFitParms::fillParams(){
  
//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   StringPairMap::const_iterator stringPairIter;

//  1. fill magnitudes and phases 
   for (int ui=paramEnumChic1ToKpKmPiGam::ChiGam; ui<paramEnumChic1ToKpKmPiGam::nAmps; ui++){

     std::string theAmpString=paramEnumChic1ToKpKmPiGam::name(ui);
     std::string theMagSuffix=theAmpString+"Mag";
     std::vector< boost::shared_ptr<const JPCLS> >  theJPCLSs=_fitParamsKpKmPiGam.jpclsVec(ui);

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagVal;
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErr;
     fillAmps(theJPCLSs, theMagSuffix, currentMagVal, currentMagErr);

     if (currentMagVal.size()>0 && currentMagErr.size()>0){
       _paramVal.Mags[ui]=currentMagVal;
       _paramErr.Mags[ui]=currentMagErr;      
     }

     std::string thePhiSuffix=theAmpString+"Phi";
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiVal;
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErr;
     fillAmps(theJPCLSs, thePhiSuffix, currentPhiVal, currentPhiErr);

     if (currentPhiVal.size()>0 && currentPhiErr.size()>0){
       _paramVal.Phis[ui]=currentPhiVal;
       _paramErr.Phis[ui]=currentPhiErr;      
     }
   }

//  2. fill masses and width 
   for (int ui=paramEnumChic1ToKpKmPiGam::nAmps; ui<paramEnumChic1ToKpKmPiGam::nMasses; ui++){
     std::string massName=paramEnumChic1ToKpKmPiGam::name(ui);
     std::string theMassSuffix=massName+"Mass";

     fillParameter(_paramVal.Masses, _paramErr.Masses, theMassSuffix, ui);

     std::string theWidthSuffix=massName+"Width";
     fillParameter(_paramVal.Widths, _paramErr.Widths, theWidthSuffix, ui);
   }

//  3. fill gFactors
   for (int ui=paramEnumChic1ToKpKmPiGam::nMasses; ui<paramEnumChic1ToKpKmPiGam::ngFactors; ui++){
     std::string gFactorName=paramEnumChic1ToKpKmPiGam::name(ui);
     std::string thegFactorSuffix=gFactorName+"gFactor";

     fillParameter(_paramVal.gFactors, _paramErr.gFactors, thegFactorSuffix, ui);

   }

//  3. fill other parameter
   for (int ui=paramEnumChic1ToKpKmPiGam::ngFactors; ui<paramEnumChic1ToKpKmPiGam::nOthers; ui++){
     std::string otherName=paramEnumChic1ToKpKmPiGam::name(ui);
     std::string theOtherSuffix=otherName+"Other";

     fillParameter(_paramVal.otherParams, _paramErr.otherParams, theOtherSuffix, ui);
   } 
}



void StreamChic1ToKpKmPiGamFitParms::fillAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& valMap , std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMap){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;

  for ( itJPCLS=theJPCLSs.begin(); itJPCLS!=theJPCLSs.end(); ++itJPCLS){

    std::string theKey=(*itJPCLS)->name()+suffix;  
    stringPairIter=_stringPairMap.find(theKey);
    
    if ( stringPairIter != _stringPairMap.end() ){
      double val=stringPairIter->second.first;
      double err=stringPairIter->second.second;
      valMap[(*itJPCLS)]=val;
      errMap[(*itJPCLS)]=err;
    } 

  }
}


void StreamChic1ToKpKmPiGamFitParms::fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index){

  StringPairMap::const_iterator stringPairIter=_stringPairMap.find(suffix);
  if ( stringPairIter != _stringPairMap.end() ){
    theValMap[index]=stringPairIter->second.first;
    theErrMap[index]=stringPairIter->second.second;    
  }

}
