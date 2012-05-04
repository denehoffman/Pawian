#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FitParamsBase.hh"
//#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParamsBase::FitParamsBase()
{
}

FitParamsBase::FitParamsBase(fitParams& theStartparams, fitParams& theErrorparams) :
  _startParams(theStartparams),
  _errorParams(theErrorparams)
{
  // fill  _ampIdx
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& startMagMap=_startParams.Mags;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){
    _ampIdx.push_back(itMagMap->first);
  }

  // fill  _ampLamLamIdx
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& startLamLamMagMap=_startParams.MagLamLams;
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::iterator itLamLamMagMap;
  for (itLamLamMagMap=startLamLamMagMap.begin(); itLamLamMagMap!=startLamLamMagMap.end(); ++itLamLamMagMap){
    _ampLamLamIdx.push_back(itLamLamMagMap->first);
  }


  // fill _massIdx
  std::map<int, double> startMassMap=_startParams.Masses; 
  std::map<int, double>::iterator itMass;
  for (itMass=startMassMap.begin(); itMass!=startMassMap.end(); ++itMass){
    _massIdx.push_back(itMass->first);
  }
  
  // fill _widthIdx
  std::map<int, double> startWidthMap=_startParams.Widths; 
  std::map<int, double>::iterator itWidth;
  for (itWidth=startWidthMap.begin(); itWidth!=startWidthMap.end(); ++itWidth){
    _widthIdx.push_back(itWidth->first);
  }

  // fill _gFactorIdx
  std::map<int, double> startgFactorMap=_startParams.gFactors; 
  std::map<int, double>::iterator itgFactor;
  for (itgFactor=startgFactorMap.begin(); itgFactor!=startgFactorMap.end(); ++itgFactor){
    _gFactorIdx.push_back(itgFactor->first);
  }

  // fill _otherIdx
  std::map<int, double> startOtherMap=_startParams.otherParams; 
  std::map<int, double>::iterator itOther;
  for (itOther=startOtherMap.begin(); itOther!=startOtherMap.end(); ++itOther){
    _otherIdx.push_back(itOther->first);
  }
}

FitParamsBase::~FitParamsBase()
{
}

std::vector< boost::shared_ptr<const JPCLS> > FitParamsBase::jpclsVec(int index){
  int nOfAmps=int(_jpclsMap.size());
  if ( index > nOfAmps ){
    Alert << "index " << index 
	  << "  number of amplitude parametes paramEnumKpKmPiGam::numAmps = " 
	  << nOfAmps << endmsg;
    exit(0);
  }

  return _jpclsMap[index]; 

}

std::vector< boost::shared_ptr<const JPClamlam> >  FitParamsBase::jpcLamLamVec(int index){
  int nOfAmps=int(_jpcLamLamMap.size());
  if ( index > nOfAmps ){
    Alert << "index " << index 
	  << "  number of amplitude parametes paramEnumKpKmPiGam::numAmps = " 
	  << nOfAmps << endmsg;
    exit(0);
  }

  return _jpcLamLamMap[index];
}


int FitParamsBase::setAmpParams( std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& ampMap, const std::vector<double>& par, int counter, int index){

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > currentStates= _jpclsMap[index];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMap;

  for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
    double currentPar=par[resultCount];
    resultCount++;
    currentMap[(*itJPCLS)]=currentPar;
  }

  ampMap[index]=currentMap;  
  return resultCount;
}

int FitParamsBase::setAmpLamLamParams( std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& ampLamLamMap, const std::vector<double>& par, int counter, int index){

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPClamlam> >::const_iterator itJPClamlam;
  std::vector< boost::shared_ptr<const JPClamlam> > currentStates= _jpcLamLamMap[index];
  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMap;

  for ( itJPClamlam=currentStates.begin(); itJPClamlam!=currentStates.end(); ++itJPClamlam){
    double currentPar=par[resultCount];
    resultCount++;
    currentMap[(*itJPClamlam)]=currentPar;
  }

  ampLamLamMap[index]=currentMap;  
  return resultCount;
}



int FitParamsBase::setSingleParams(std::map<int, double>& theMap, const std::vector<double>& par, int counter, int index){
  int resultCount=counter;
  theMap[index]=par[resultCount];
  resultCount++;
  return resultCount;
}

void FitParamsBase::setMnUsrParams(MnUserParameters& upar){

  // 1.: set magnitudes of all lamlam amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& startLamLamMagMap=_startParams.MagLamLams;
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& errLamLamMagMap=_errorParams.PhiLamLams;

  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::iterator itLamLamMagMap;
  for (itLamLamMagMap=startLamLamMagMap.begin(); itLamLamMagMap!=startLamLamMagMap.end(); ++itLamLamMagMap){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& errLamLamMags=errLamLamMagMap[itLamLamMagMap->first];
    
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLamMag;
    for (itLamLamMag=itLamLamMagMap->second.begin(); itLamLamMag!=itLamLamMagMap->second.end(); ++itLamLamMag){
      int theIndex=itLamLamMagMap->first;
      boost::shared_ptr<const JPClamlam> theJPCLamLam=itLamLamMag->first;
      double theStartVal=itLamLamMag->second;
      double theErrVal=errLamLamMags[theJPCLamLam];
      
      //now fill the fitParameterMap
      std::string magStr=theJPCLamLam->name()+ampName(theIndex)+"Mag";
      
      double magMin=theStartVal-3.*theErrVal;
      if (magMin<0.) magMin=0.;
      
      upar.Add(magStr, theStartVal, theErrVal, magMin, theStartVal+3.*theErrVal);
    }
    
  }

   // 2.: set phases of all lam lam amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& startLamLamPhiMap=_startParams.PhiLamLams;
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& errLamLamPhiMap=_errorParams.PhiLamLams;

  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::iterator itLamLamPhiMap;
  for (itLamLamPhiMap=startLamLamPhiMap.begin(); itLamLamPhiMap!=startLamLamPhiMap.end(); ++itLamLamPhiMap){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& errLamLamPhis=errLamLamPhiMap[itLamLamPhiMap->first];
    
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLamPhi;
    for (itLamLamPhi=itLamLamPhiMap->second.begin(); itLamLamPhi!=itLamLamPhiMap->second.end(); ++itLamLamPhi){
      int theIndex=itLamLamPhiMap->first;
      boost::shared_ptr<const JPClamlam> theJPCLamLam=itLamLamPhi->first;
      double theStartVal=itLamLamPhi->second;
      double theErrVal=errLamLamPhis[theJPCLamLam];
      
      //now fill the fitParameterMap
      std::string phiStr=theJPCLamLam->name()+ampName(theIndex)+"Phi";
      
      double phiMin=theStartVal-3.*theErrVal;
      if (phiMin<0.) phiMin=0.;
      
      upar.Add(phiStr, theStartVal, theErrVal, -3.*M_PI, 3.*M_PI);
    }
    
  }


  // 3.: set magnitudes of all amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& startMagMap=_startParams.Mags;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& errMagMap=_errorParams.Mags;

  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMags=errMagMap[itMagMap->first];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagMap->second.begin(); itMag!=itMagMap->second.end(); ++itMag){
      int theIndex=itMagMap->first;
      boost::shared_ptr<const JPCLS> theJPCLS=itMag->first;
      double theStartVal=itMag->second;
      double theErrVal=errMags[theJPCLS];
      
      //now fill the fitParameterMap
      std::string magStr=theJPCLS->name()+ampName(theIndex)+"Mag";
      
      double magMin=theStartVal-3.*theErrVal;
      if (magMin<0.) magMin=0.;
      
      upar.Add(magStr, theStartVal, theErrVal, magMin, theStartVal+3.*theErrVal);
    }
    
  }

   // 4.: set phases of all amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& startPhiMap=_startParams.Phis;
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& errPhiMap=_errorParams.Phis;

  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itPhiMap;
  for (itPhiMap=startPhiMap.begin(); itPhiMap!=startPhiMap.end(); ++itPhiMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errPhis=errPhiMap[itPhiMap->first];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPhi;
    for (itPhi=itPhiMap->second.begin(); itPhi!=itPhiMap->second.end(); ++itPhi){
      int theIndex=itPhiMap->first;
      boost::shared_ptr<const JPCLS> theJPCLS=itPhi->first;
      double theStartVal=itPhi->second;
      double theErrVal=errPhis[theJPCLS];
      
      //now fill the fitParameterMap
      std::string phiStr=theJPCLS->name()+ampName(theIndex)+"Phi";
      
      double phiMin=theStartVal-3.*theErrVal;
      if (phiMin<0.) phiMin=0.;
      
      upar.Add(phiStr, theStartVal, theErrVal, -3.*M_PI, 3.*M_PI);
    }
    
  }

   // 5.: set all masses
  std::map<int, double>& startMassMap= _startParams.Masses;
  std::map<int, double>& errMassMap= _errorParams.Masses;

  std::map<int, double>::iterator itMass;
  for (itMass=startMassMap.begin(); itMass!=startMassMap.end(); ++itMass){
     int theIndex=itMass->first;
     double theStartVal=itMass->second;
     double theErrVal=errMassMap[theIndex];
 
      std::string massStr=massName(theIndex)+"Mass";    

      double massMin=theStartVal-5.*theErrVal;
      if (massMin<0.) massMin=0.;

      double massMax=theStartVal+5.*theErrVal;
      upar.Add(massStr, theStartVal, theErrVal, massMin, massMax);
  }

   // 6.: set all widths
  std::map<int, double>& startWidthMap= _startParams.Widths;
  std::map<int, double>& errWidthMap= _errorParams.Widths;

  std::map<int, double>::iterator itWidth;
  for (itWidth=startWidthMap.begin(); itWidth!=startWidthMap.end(); ++itWidth){
     int theIndex=itWidth->first;
     double theStartVal=itWidth->second;
     double theErrVal=errWidthMap[theIndex];
 
      std::string widthStr=widthName(theIndex)+"Width";    

      double widthMin=0.;
      double widthMax=2.*theStartVal;
      upar.Add(widthStr, theStartVal, theErrVal, widthMin, widthMax);
  } 

   // 7.: set all gFactors
  std::map<int, double>& startgFactorMap= _startParams.gFactors;
  std::map<int, double>& errgFactorMap= _errorParams.gFactors;

  std::map<int, double>::iterator itgFactor;
  for (itgFactor=startgFactorMap.begin(); itgFactor!=startgFactorMap.end(); ++itgFactor){
     int theIndex=itgFactor->first;
     double theStartVal=itgFactor->second;
     double theErrVal=errgFactorMap[theIndex];
 
      std::string gFactorStr=gFactorName(theIndex)+"gFactor";    

      double gFactorMin=theStartVal-5.*theErrVal;
      if (gFactorMin<0.) gFactorMin=0.;
      double gFactorMax=theStartVal+5.*theErrVal;
      upar.Add(gFactorStr, theStartVal, theErrVal, gFactorMin, gFactorMax);
  } 

   // 8.: set all other parameters
  std::map<int, double>& startotherParamsMap= _startParams.otherParams;
  std::map<int, double>& errotherParamsMap= _errorParams.otherParams;

  std::map<int, double>::iterator itotherParams;
  for (itotherParams=startotherParamsMap.begin(); itotherParams!=startotherParamsMap.end(); ++itotherParams){
     int theIndex=itotherParams->first;
     double theStartVal=itotherParams->second;
     double theErrVal=errotherParamsMap[theIndex];
 
      std::string otherParamsStr=otherName(theIndex)+"OtherParams";    

      double otherParamsMin=theStartVal-5.*theErrVal;
      if (otherParamsMin<0.) otherParamsMin=0.;
      double otherParamsMax=theStartVal+5.*theErrVal;
      upar.Add(otherParamsStr, theStartVal, theErrVal, otherParamsMin, otherParamsMax);
  } 
}



void FitParamsBase::printParams(fitParams& theParams){

  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::const_iterator itLamLamAmps;

  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLamLamMap;

   for (itLamLamAmps=theParams.MagLamLams.begin(); itLamLamAmps!=theParams.MagLamLams.end(); ++itLamLamAmps){
     std::string currentAmpName=ampName(itLamLamAmps->first);

     std::string currentMagName=currentAmpName+"Mag";
     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=itLamLamAmps->second;    
     for ( itJPCLamLamMap=currentAmp.begin(); itJPCLamLamMap!=currentAmp.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentMagName;
        Info<< currentName << "\t" << itJPCLamLamMap->second << endmsg;
     }
   }

  for (itLamLamAmps=theParams.PhiLamLams.begin(); itLamLamAmps!=theParams.PhiLamLams.end(); ++itLamLamAmps){

    std::string currentAmpName=ampName(itLamLamAmps->first);
    std::string currentPhiName=currentAmpName+"Phi";

    for (itJPCLamLamMap=itLamLamAmps->second.begin(); itJPCLamLamMap!=itLamLamAmps->second.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentPhiName;
      Info<< currentName <<"\t" << itJPCLamLamMap->second << endmsg;
    }
  }


  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::const_iterator itAmps;

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLSMap;

   for (itAmps=theParams.Mags.begin(); itAmps!=theParams.Mags.end(); ++itAmps){
     std::string currentAmpName=ampName(itAmps->first);

     std::string currentMagName=currentAmpName+"Mag";
     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=itAmps->second;    
     for ( itJPCLSMap=currentAmp.begin(); itJPCLSMap!=currentAmp.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+currentMagName;
        Info<< currentName << "\t" << itJPCLSMap->second << endmsg;
     }
   }

  for (itAmps=theParams.Phis.begin(); itAmps!=theParams.Phis.end(); ++itAmps){

    std::string currentAmpName=ampName(itAmps->first);
    std::string currentPhiName=currentAmpName+"Phi";

    for (itJPCLSMap=itAmps->second.begin(); itJPCLSMap!=itAmps->second.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+currentPhiName;
      Info<< currentName <<"\t" << itJPCLSMap->second << endmsg;
    }
  }

  std::map<int, double>::const_iterator itSinglePar;

  if (theParams.Masses.size()>0){
    for (itSinglePar=theParams.Masses.begin(); itSinglePar!=theParams.Masses.end(); ++itSinglePar){
      Info << massName(itSinglePar->first) << "Mass\t" << itSinglePar->second;
    }
  }

  if (theParams.Widths.size()>0){
    for (itSinglePar=theParams.Widths.begin(); itSinglePar!=theParams.Widths.end(); ++itSinglePar){
      Info << widthName(itSinglePar->first) << "Width\t" << itSinglePar->second;
    }
  }

  if (theParams.gFactors.size()>0){
    for (itSinglePar=theParams.gFactors.begin(); itSinglePar!=theParams.gFactors.end(); ++itSinglePar){
      Info << gFactorName(itSinglePar->first) << "gFactor\t" << itSinglePar->second;
    }
  }

  if (theParams.otherParams.size()>0){
    for (itSinglePar=theParams.otherParams.begin(); itSinglePar!=theParams.otherParams.end(); ++itSinglePar){
      Info << otherName(itSinglePar->first) << "\t" << itSinglePar->second;
    }
  }

}

void FitParamsBase::dumpParams(std::ostream& os, fitParams& theVals,  fitParams& theErrs){

  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::const_iterator itLamLamAmps;

  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLamLamMap;

   for (itLamLamAmps=theVals.MagLamLams.begin(); itLamLamAmps!=theVals.MagLamLams.end(); ++itLamLamAmps){
     std::string currentAmpName=ampName(itLamLamAmps->first);
     std::string currentMagName=currentAmpName+"Mag";

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentErrMap=theErrs.MagLamLams[itLamLamAmps->first];

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=itLamLamAmps->second;    
     for ( itJPCLamLamMap=currentAmp.begin(); itJPCLamLamMap!=currentAmp.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentMagName;
       
       os << currentName << "\t" << itJPCLamLamMap->second << "\t" << currentErrMap[itJPCLamLamMap->first] << std::endl;
     }
   }

   for (itLamLamAmps=theVals.PhiLamLams.begin(); itLamLamAmps!=theVals.PhiLamLams.end(); ++itLamLamAmps){
     std::string currentAmpName=ampName(itLamLamAmps->first);
     std::string currentMagName=currentAmpName+"Phi";

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentErrMap=theErrs.PhiLamLams[itLamLamAmps->first];

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=itLamLamAmps->second;    
     for ( itJPCLamLamMap=currentAmp.begin(); itJPCLamLamMap!=currentAmp.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentMagName;
       
       os << currentName << "\t" << itJPCLamLamMap->second << "\t" << currentErrMap[itJPCLamLamMap->first] << std::endl;
     }
   }



  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::const_iterator itAmps;

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLSMap;

   for (itAmps=theVals.Mags.begin(); itAmps!=theVals.Mags.end(); ++itAmps){
     std::string currentAmpName=ampName(itAmps->first);
     std::string currentMagName=currentAmpName+"Mag";

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentErrMap=theErrs.Mags[itAmps->first];

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=itAmps->second;    
     for ( itJPCLSMap=currentAmp.begin(); itJPCLSMap!=currentAmp.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+currentMagName;
       
       os << currentName << "\t" << itJPCLSMap->second << "\t" << currentErrMap[itJPCLSMap->first] << std::endl;
     }
   }

   for (itAmps=theVals.Phis.begin(); itAmps!=theVals.Phis.end(); ++itAmps){
     std::string currentAmpName=ampName(itAmps->first);
     std::string currentMagName=currentAmpName+"Phi";

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentErrMap=theErrs.Phis[itAmps->first];

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=itAmps->second;    
     for ( itJPCLSMap=currentAmp.begin(); itJPCLSMap!=currentAmp.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+currentMagName;
       
       os << currentName << "\t" << itJPCLSMap->second << "\t" << currentErrMap[itJPCLSMap->first] << std::endl;
     }
   }

  std::map<int, double>::const_iterator itSinglePar;
  if (theVals.Masses.size()>0 && theErrs.Masses.size()==theVals.Masses.size()){
    for (itSinglePar=theVals.Masses.begin(); itSinglePar!=theVals.Masses.end(); ++itSinglePar){
      std::string currentName=massName(itSinglePar->first)+"Mass";
      os << currentName << "\t" << itSinglePar->second << "\t" << theErrs.Masses[itSinglePar->first] << std::endl;
    }
  }

  if (theVals.Widths.size()>0 && theErrs.Widths.size()==theVals.Widths.size()){
    for (itSinglePar=theVals.Widths.begin(); itSinglePar!=theVals.Widths.end(); ++itSinglePar){
      std::string currentName=widthName(itSinglePar->first)+"Width";
      os << currentName << "\t" << itSinglePar->second << "\t" << theErrs.Widths[itSinglePar->first] << std::endl;
    }
  }

  if (theVals.gFactors.size()>0 && theErrs.gFactors.size()==theVals.gFactors.size()){
    for (itSinglePar=theVals.gFactors.begin(); itSinglePar!=theVals.gFactors.end(); ++itSinglePar){
      std::string currentName=gFactorName(itSinglePar->first)+"gFactor";
      os << currentName << "\t" << itSinglePar->second << "\t" << theErrs.gFactors[itSinglePar->first] << std::endl;
    }
  }

  if (theVals.otherParams.size()>0 && theErrs.otherParams.size()==theVals.otherParams.size()){
    for (itSinglePar=theVals.otherParams.begin(); itSinglePar!=theVals.otherParams.end(); ++itSinglePar){
      std::string currentName=otherName(itSinglePar->first)+"Other";
      os << currentName << "\t" << itSinglePar->second << "\t" << theErrs.otherParams[itSinglePar->first] << std::endl;
    }
  }
}

fitParams FitParamsBase::getFitParamVal(const std::vector<double>& par){
  fitParams result=_startParams;
  int counter=0;

  // 1.: magnitudes of all lamlam amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& magLamLamMap=result.MagLamLams;

  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::iterator itMagLamLamMap;
  for (itMagLamLamMap=magLamLamMap.begin(); itMagLamLamMap!=magLamLamMap.end(); ++itMagLamLamMap){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagLamLamMap->second.begin(); itMag!=itMagLamLamMap->second.end(); ++itMag){
      itMag->second=par[counter];
      counter++;

    }
  }

  // 2.: phases of all lam lam amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& phiLamLamMap=result.PhiLamLams;
  
  std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >::iterator itPhiLamLamMap;
  for (itPhiLamLamMap=phiLamLamMap.begin(); itPhiLamLamMap!=phiLamLamMap.end(); ++itPhiLamLamMap){
    
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itPhi;
    for (itPhi=itPhiLamLamMap->second.begin(); itPhi!=itPhiLamLamMap->second.end(); ++itPhi){
      itPhi->second=par[counter];
      counter++;
    }
  }


  // 3.: magnitudes of all amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& magMap=result.Mags;

  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itMagMap;
  for (itMagMap=magMap.begin(); itMagMap!=magMap.end(); ++itMagMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagMap->second.begin(); itMag!=itMagMap->second.end(); ++itMag){
      itMag->second=par[counter];
      counter++;

    }
  }

  // 4.: phases of all amplitudes
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& phiMap=result.Phis;
  
  std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >::iterator itPhiMap;
  for (itPhiMap=phiMap.begin(); itPhiMap!=phiMap.end(); ++itPhiMap){
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPhi;
    for (itPhi=itPhiMap->second.begin(); itPhi!=itPhiMap->second.end(); ++itPhi){
      itPhi->second=par[counter];
      counter++;
    }
  }

  // 5.: masses
  std::map<int, double>& massMap= result.Masses;
  std::map<int, double>::iterator itMass;
  for (itMass=massMap.begin(); itMass!=massMap.end(); ++itMass){
      itMass->second=par[counter];
      counter++;
 }

  // 6.: widths
  std::map<int, double>& widthMap= result.Widths;
  std::map<int, double>::iterator itWidth;
  for (itWidth=widthMap.begin(); itWidth!=widthMap.end(); ++itWidth){
      itWidth->second=par[counter];
      counter++;
 }

  // 7.: gFactors
  std::map<int, double>& gFactorMap= result.gFactors;
  std::map<int, double>::iterator itgFactor;
  for (itgFactor=gFactorMap.begin(); itgFactor!=gFactorMap.end(); ++itgFactor){
      itgFactor->second=par[counter];
      counter++;
 }

  // 8.: others
  std::map<int, double>& otherMap= result.otherParams;
  std::map<int, double>::iterator itOther;
  for (itOther=otherMap.begin(); itOther!=otherMap.end(); ++itOther){
      itOther->second=par[counter];
      counter++;
 }
  return result;
}
