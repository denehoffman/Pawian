#include <fstream>

#include "PwaUtils/StreamFitParmsBase.hh"
#include "ErrLogger/ErrLogger.hh"

StreamFitParmsBase::StreamFitParmsBase(std::string& filePath, boost::shared_ptr<FitParamsBase> fitParamsBasePtr) :
  AbsFitParamStreamer(filePath),
  _fitParamsBasePtr(fitParamsBasePtr)
{
  fillParams();

}

StreamFitParmsBase::~StreamFitParmsBase(){;}

void StreamFitParmsBase::fillParams(){
  
//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   StringPairMap::const_iterator stringPairIter;

//  0. fill lam lam magnitudes and phases 
   for (int ui=_fitParamsBasePtr->ampLamLamIdxMin(); ui<=_fitParamsBasePtr->ampLamLamIdxMax(); ui++){

     std::string theAmpString=_fitParamsBasePtr->ampLamLamName(ui);
     std::string theMagSuffix=theAmpString+"Mag";
     std::vector< boost::shared_ptr<const JPClamlam> >  theJPCLamLams=_fitParamsBasePtr->jpcLamLamVec(ui);

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagVal;
     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMagErr;
     fillLamLamAmps(theJPCLamLams, theMagSuffix, currentMagVal, currentMagErr);

     if (currentMagVal.size()>0 && currentMagErr.size()>0){
       _paramVal.MagLamLams[ui]=currentMagVal;
       _paramErr.MagLamLams[ui]=currentMagErr;      
     }

     std::string thePhiSuffix=theAmpString+"Phi";
     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiVal;
     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentPhiErr;
     fillLamLamAmps(theJPCLamLams, thePhiSuffix, currentPhiVal, currentPhiErr);

     if (currentPhiVal.size()>0 && currentPhiErr.size()>0){
       _paramVal.PhiLamLams[ui]=currentPhiVal;
       _paramErr.PhiLamLams[ui]=currentPhiErr;      
     }
   }

//  1. fill magnitudes and phases 
   for (int ui=_fitParamsBasePtr->ampIdxMin(); ui<=_fitParamsBasePtr->ampIdxMax(); ui++){

     std::string theAmpString=_fitParamsBasePtr->ampName(ui);
     std::string theMagSuffix=theAmpString+"Mag";
     std::vector< boost::shared_ptr<const JPCLS> >  theJPCLSs=_fitParamsBasePtr->jpclsVec(ui);

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
   for (int ui=_fitParamsBasePtr->massIdxMin(); ui<=_fitParamsBasePtr->massIdxMax(); ui++){
     std::string massName=_fitParamsBasePtr->massName(ui);
     std::string theMassSuffix=massName+"Mass";

     fillParameter(_paramVal.Masses, _paramErr.Masses, theMassSuffix, ui);

     std::string theWidthSuffix=massName+"Width";
     fillParameter(_paramVal.Widths, _paramErr.Widths, theWidthSuffix, ui);
   }

//  3. fill gFactors
   for (int ui=_fitParamsBasePtr->gFactorIdxMin(); ui<=_fitParamsBasePtr->gFactorIdxMax(); ui++){
     std::string gFactorName=_fitParamsBasePtr->gFactorName(ui);
     std::string thegFactorSuffix=gFactorName+"gFactor";

     fillParameter(_paramVal.gFactors, _paramErr.gFactors, thegFactorSuffix, ui);

   }

//  3. fill other parameter
   for (int ui=_fitParamsBasePtr->otherIdxMin(); ui<=_fitParamsBasePtr->otherIdxMax(); ui++){
     std::string otherName=_fitParamsBasePtr->otherName(ui);
     std::string theOtherSuffix=otherName+"Other";

     fillParameter(_paramVal.otherParams, _paramErr.otherParams, theOtherSuffix, ui);
   } 
}



void StreamFitParmsBase::fillAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& valMap , std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMap){

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

void StreamFitParmsBase::fillLamLamAmps(std::vector< boost::shared_ptr<const JPClamlam> >& theJPCLamLams, std::string& suffix, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& valMap , std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& errMap){
 
  std::vector< boost::shared_ptr<const JPClamlam> >::const_iterator itJPCLamLam;
  StringPairMap::const_iterator stringPairIter;

  for ( itJPCLamLam=theJPCLamLams.begin(); itJPCLamLam!=theJPCLamLams.end(); ++itJPCLamLam){

    std::string theKey=(*itJPCLamLam)->name()+suffix;  
    stringPairIter=_stringPairMap.find(theKey);
    
    if ( stringPairIter != _stringPairMap.end() ){
      double val=stringPairIter->second.first;
      double err=stringPairIter->second.second;
      valMap[(*itJPCLamLam)]=val;
      errMap[(*itJPCLamLam)]=err;
    } 

  }
}

void StreamFitParmsBase::fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index){

  StringPairMap::const_iterator stringPairIter=_stringPairMap.find(suffix);
  if ( stringPairIter != _stringPairMap.end() ){
    theValMap[index]=stringPairIter->second.first;
    theErrMap[index]=stringPairIter->second.second;    
  }

}
