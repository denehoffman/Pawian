#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FitParamsBase.hh"
//#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParamsBase::FitParamsBase() :
  _magSuffix("Mag")
  ,_phiSuffix("Phi")
  ,_massSuffix("Mass")
  ,_widthSuffix("Width")
  ,_gFactorSuffix("gFactor")
  ,_otherParamSuffix("Other")
{
}


FitParamsBase::~FitParamsBase()
{
}

// int FitParamsBase::setAmpParams( mapStrJPCLS& ampMap, const std::vector<double>& par, int counter, std::string& key){

//   int resultCount=counter;
//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMap=
//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   std::vector< boost::shared_ptr<const JPCLS> > currentStates= _jpclsMap[index];

//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMap;

//   for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
//     double currentPar=par[resultCount];
//     resultCount++;
//     currentMap[(*itJPCLS)]=currentPar;
//   }

//   ampMap[key]=currentMap;  
//   return resultCount;
// }

// int FitParamsBase::setAmpLamLamParams( std::map<int, std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > >& ampLamLamMap, const std::vector<double>& par, int counter, int index){

//   int resultCount=counter;
//   std::vector< boost::shared_ptr<const JPClamlam> >::const_iterator itJPClamlam;
//   std::vector< boost::shared_ptr<const JPClamlam> > currentStates= _jpcLamLamMap[index];
//   std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentMap;

//   for ( itJPClamlam=currentStates.begin(); itJPClamlam!=currentStates.end(); ++itJPClamlam){
//     double currentPar=par[resultCount];
//     resultCount++;
//     currentMap[(*itJPClamlam)]=currentPar;
//   }

//   ampLamLamMap[index]=currentMap;  
//   return resultCount;
// }



// int FitParamsBase::setSingleParams(std::map<int, double>& theMap, const std::vector<double>& par, int counter, int index){
//   int resultCount=counter;
//   theMap[index]=par[resultCount];
//   resultCount++;
//   return resultCount;
// }

void FitParamsBase::setMnUsrParams(MnUserParameters& upar, fitParams& theValParams, fitParams& theErrParams){

  // 1.: set magnitudes of all lamlam amplitudes
  setMnUsrParamsJPCLamLam(upar, theValParams.MagLamLams, theErrParams.MagLamLams, _magSuffix);

  // 2.: set phases of all lam lam amplitudes
  setMnUsrParamsJPCLamLam(upar, theValParams.PhiLamLams, theErrParams.PhiLamLams, _phiSuffix);

  // 3.: set magnitudes of all amplitudes
  setMnUsrParamsJPCLS(upar, theValParams.Mags, theErrParams.Mags, _magSuffix);

   // 4.: set phases of all amplitudes
  setMnUsrParamsJPCLS(upar, theValParams.Phis, theErrParams.Phis, _phiSuffix);


   // 5.: set all masses
  setMnUsrParamsDouble(upar, theValParams.Masses, theErrParams.Masses, _massSuffix);

   // 6.: set all widths
   setMnUsrParamsDouble(upar, theValParams.Widths, theErrParams.Widths, _widthSuffix);
 
   // 7.: set all gFactors 
  setMnUsrParamsDouble(upar, theValParams.gFactors, theErrParams.gFactors, _gFactorSuffix);

   // 8.: set all other parameters
  setMnUsrParamsDouble(upar, theValParams.otherParams, theErrParams.otherParams, _otherParamSuffix);

}



void FitParamsBase::printParams(fitParams& theParams){

  printJPCLamLamParams(theParams, theParams.MagLamLams, _magSuffix);
  printJPCLamLamParams(theParams, theParams.PhiLamLams, _phiSuffix);

  printJPCLSParams(theParams, theParams.Mags, _magSuffix);
  printJPCLSParams(theParams, theParams.Phis, _phiSuffix);

  printDoubleParams(theParams, theParams.Masses, _massSuffix);
  printDoubleParams(theParams, theParams.Widths, _widthSuffix);
  printDoubleParams(theParams, theParams.gFactors, _gFactorSuffix);
  printDoubleParams(theParams, theParams.otherParams, _otherParamSuffix);
}

void FitParamsBase::dumpParams(std::ostream& os, fitParams& theVals,  fitParams& theErrs){

  dumpLamLamParams(os, theVals.MagLamLams, theErrs.MagLamLams, _magSuffix);
  dumpLamLamParams(os, theVals.PhiLamLams, theErrs.PhiLamLams, _phiSuffix);

  dumpJPCLSParams(os, theVals.Mags, theErrs.Mags,  _magSuffix);
  dumpJPCLSParams(os, theVals.Phis, theErrs.Phis,  _phiSuffix);

  dumpDoubleParams(os, theVals.Masses, theErrs.Masses, _massSuffix);
  dumpDoubleParams(os, theVals.Widths, theErrs.Widths, _widthSuffix);
  dumpDoubleParams(os, theVals.gFactors, theErrs.gFactors, _gFactorSuffix);
  dumpDoubleParams(os, theVals.otherParams, theErrs.otherParams, _otherParamSuffix);

}

void FitParamsBase::getFitParamVal(const std::vector<double>& par, fitParams& theParams){

  int counter=0;
  getFitParamValJPCLamLam(par, theParams.MagLamLams, counter);
  getFitParamValJPCLamLam(par, theParams.PhiLamLams, counter);
  getFitParamValJPCLS(par, theParams.Mags, counter);
  getFitParamValJPCLS(par, theParams.Phis, counter);
  getFitParamValDouble(par, theParams.Masses, counter);
  getFitParamValDouble(par, theParams.Widths, counter);
  getFitParamValDouble(par, theParams.gFactors, counter);
  getFitParamValDouble(par, theParams.otherParams, counter);

}







void FitParamsBase::setMnUsrParamsJPCLamLam(MnUserParameters& upar, mapStrJPCLamLam& startLamLamMagMap, mapStrJPCLamLam& errLamLamMagMap, const std::string& suffix){
 
  mapStrJPCLamLam::iterator itLamLamMagMap;
  for (itLamLamMagMap=startLamLamMagMap.begin(); itLamLamMagMap!=startLamLamMagMap.end(); ++itLamLamMagMap){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& errLamLamMags= errLamLamMagMap[itLamLamMagMap->first];
    
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLamMag;
    for (itLamLamMag=itLamLamMagMap->second.begin(); itLamLamMag!=itLamLamMagMap->second.end(); ++itLamLamMag){
      
      boost::shared_ptr<const JPClamlam> theJPCLamLam=itLamLamMag->first;
      double theStartVal=itLamLamMag->second;
      double theErrVal=errLamLamMags[theJPCLamLam];
      
      //now fill the fitParameterMap
      std::string magStr=theJPCLamLam->name()+itLamLamMagMap->first+suffix;

      double valMin=0.;
      double valMax=0.;
      if (suffix==_phiSuffix){
	valMin=-4.*M_PI;
	valMax=4.*M_PI;
      }
      else{
	valMin=theStartVal-3.*theErrVal;
	if (valMin<0.) valMin=0.;
	valMax=theStartVal+3.*theErrVal;
      }
      
      upar.Add(magStr, theStartVal, theErrVal, valMin, valMax);
    }
    
  }

}

void FitParamsBase::setMnUsrParamsJPCLS(MnUserParameters& upar, mapStrJPCLS& startMagMap, mapStrJPCLS& errMagMap, const std::string& suffix){

  mapStrJPCLS::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMags=errMagMap[itMagMap->first];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagMap->second.begin(); itMag!=itMagMap->second.end(); ++itMag){

      boost::shared_ptr<const JPCLS> theJPCLS=itMag->first;
      double theStartVal=itMag->second;
      double theErrVal=errMags[theJPCLS];
      if (theErrVal<=0.) theErrVal=0.1;
      
      //now fill the fitParameterMap
      std::string magStr=theJPCLS->name()+itMagMap->first+suffix;

      double valMin=0.;
      double valMax=0.;
      if (suffix==_phiSuffix){
	valMin=-4.*M_PI;
	valMax=4.*M_PI;
      }
      else{
	valMin=theStartVal-3.*theErrVal;
	if (valMin<0.) valMin=0.;
	valMax=theStartVal+3.*theErrVal;
      }      

      upar.Add(magStr, theStartVal, theErrVal, valMin, valMax);
    }
    
  }

}

void FitParamsBase::setMnUsrParamsDouble(MnUserParameters& upar, mapStrDouble& startDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix){

  mapStrDouble::iterator it;
  for (it=startDoubleMap.begin(); it!=startDoubleMap.end(); ++it){
    std::string theName =it->first + suffix;
    double theStartVal=it->second;
    double theErrVal=errDoubleMap[it->first];
    if (theErrVal<=0.) theErrVal=0.1;
    
    double minVal=theStartVal-5.*theErrVal;
    if (minVal<0.) minVal=0.;
    
    double maxVal=theStartVal+5.*theErrVal;

    // for complex fit parameter; phi component; quick workaround
    if(theName.size()>9){
      if (theName.compare(theName.size()-8, theName.size(), "PhiOther")==0){
	minVal=-4.*M_PI;
	maxVal=4.*M_PI;	    
      }
    }

    // for parameter where pos and neg values are allowed
    if(theName.size()>12){
      if (theName.compare(theName.size()-11, theName.size(), "PosNegOther")==0){
	minVal = -fabs(theStartVal)-5.*theErrVal;
	maxVal = fabs(theStartVal)+5.*theErrVal;;	    
      }
    }

   
    upar.Add(theName, theStartVal, theErrVal, minVal, maxVal);
  }

}

void FitParamsBase::printJPCLamLamParams(fitParams& theParams, mapStrJPCLamLam& lamLamMap, const std::string& suffix){

  mapStrJPCLamLam::const_iterator it;

  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLamLamMap;

   for (it=lamLamMap.begin(); it!=lamLamMap.end(); ++it){
     std::string currentAmpName=it->first;

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;    
     for ( itJPCLamLamMap=currentAmp.begin(); itJPCLamLamMap!=currentAmp.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentAmpName+suffix;
       
       Info << currentName << "\t" << itJPCLamLamMap->second << endmsg;
     }
   }
}


void FitParamsBase::printJPCLSParams(fitParams& theParams, mapStrJPCLS& JPCLSMagMap, const std::string& suffix){

  mapStrJPCLS::const_iterator it;

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLSMap;

   for (it=JPCLSMagMap.begin(); it!=JPCLSMagMap.end(); ++it){
 
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;    
     for ( itJPCLSMap=currentAmp.begin(); itJPCLSMap!=currentAmp.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+it->first+suffix;
        Info<< currentName << "\t" << itJPCLSMap->second << endmsg;
     }
   }

}

void FitParamsBase::printDoubleParams(fitParams& theParams, mapStrDouble& doubleMap, const std::string& suffix){

  mapStrDouble::const_iterator it;
  
  for (it=doubleMap.begin(); it!=doubleMap.end(); ++it){
    std::string currentName=it->first+suffix;
    Info << currentName << "\t" << it->second;
  }
  

}

void FitParamsBase::getFitParamValJPCLamLam(const std::vector<double>& par, mapStrJPCLamLam& lamLamMagMap, int& counter){
//   Info << "getFitParamValJPCLamLam par[" << counter << "]=\t" << par[counter] << endmsg;
  mapStrJPCLamLam::iterator it;
  for (it=lamLamMagMap.begin(); it!=lamLamMagMap.end(); ++it){

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itJPClamlam;
    for (itJPClamlam=it->second.begin(); itJPClamlam!=it->second.end(); ++itJPClamlam){
      itJPClamlam->second=par[counter];
      counter++;

    }
  }
}

void FitParamsBase::getFitParamValJPCLS(const std::vector<double>& par, mapStrJPCLS& jpclsMap, int& counter){
//   Info << "getFitParamValJPCLS par[" << counter << "]=\t" << par[counter] << endmsg;
  mapStrJPCLS::iterator it;
  for (it=jpclsMap.begin(); it!=jpclsMap.end(); ++it){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itjpcls;
    for (itjpcls=it->second.begin(); itjpcls!=it->second.end(); ++itjpcls){
      itjpcls->second=par[counter];
      counter++;

    }
  }
}


void FitParamsBase::getFitParamValDouble(const std::vector<double>& par, mapStrDouble& doubleMap, int& counter){
  //   Info << "getFitParamValDouble par[" << counter << "]=\t" << par[counter] << endmsg;
  mapStrDouble::iterator it;
  for (it=doubleMap.begin(); it!=doubleMap.end(); ++it){
    it->second=par[counter];
    counter++;
  }
}


void  FitParamsBase::dumpLamLamParams(std::ostream& os, mapStrJPCLamLam& valLamLamMap, mapStrJPCLamLam& errLamLamMap, const std::string& suffix){

  mapStrJPCLamLam::const_iterator it;

   for (it=valLamLamMap.begin(); it!=valLamLamMap.end(); ++it){

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentErrMap=errLamLamMap[it->first];

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;

     std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itLamLam;   
     for ( itLamLam=currentAmp.begin(); itLamLam!=currentAmp.end(); ++itLamLam){
       std::string currentName=itLamLam->first->name()+it->first+suffix;
       
       os << currentName << "\t" << itLamLam->second << "\t" << currentErrMap[itLamLam->first] << std::endl;
     }
   }

}

void FitParamsBase::dumpJPCLSParams(std::ostream& os, mapStrJPCLS& valLSMap, mapStrJPCLS& errLSMap, const std::string& suffix){

  mapStrJPCLS::const_iterator it;

   for (it=valLSMap.begin(); it!=valLSMap.end(); ++it){

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentErrMap=errLSMap[it->first];

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;

     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itLS;   
     for ( itLS=currentAmp.begin(); itLS!=currentAmp.end(); ++itLS){
       std::string currentName=itLS->first->name()+it->first+suffix;
       
       os << currentName << "\t" << itLS->second << "\t" << currentErrMap[itLS->first] << std::endl;
     }
   }
}

void FitParamsBase::dumpDoubleParams(std::ostream& os, mapStrDouble& valDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix){

  mapStrDouble::const_iterator it;

  for (it=valDoubleMap.begin(); it!=valDoubleMap.end(); ++it){
    std::string currentName=it->first+suffix;
      os << currentName << "\t" << it->second << "\t" << errDoubleMap[it->first] << std::endl;
    }

}
