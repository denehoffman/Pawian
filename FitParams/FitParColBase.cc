//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

#include <getopt.h>
#include <fstream>
#include <string>
#include <iomanip>

#include "FitParams/FitParColBase.hh"
//#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParColBase::FitParColBase() :
  _magSuffix("Mag")
  ,_phiSuffix("Phi")
  ,_massSuffix("Mass")
  ,_widthSuffix("Width")
  ,_gFactorSuffix("gFactor")
  ,_otherParamSuffix("Other")
{
}


FitParColBase::~FitParColBase()
{
}

void FitParColBase::setAbsPawianParams(std::shared_ptr<AbsPawianParameters> upar, fitParCol& theValParams, fitParCol& theErrParams){

  // 0.: set magnitudes of all JPC amplitudes
   setAbsPawianParamsJPC(upar, theValParams.MagsJPC, theErrParams.MagsJPC, _magSuffix);

  // 0a.: set phi of all JPC amplitudes
   setAbsPawianParamsJPC(upar, theValParams.PhisJPC, theErrParams.PhisJPC, _phiSuffix);

  // 1.: set magnitudes of all JPCLS amplitudes
  setAbsPawianParamsJPCLS(upar, theValParams.Mags, theErrParams.Mags, _magSuffix);

   // 2.: set phases of all JPCLS amplitudes
  setAbsPawianParamsJPCLS(upar, theValParams.Phis, theErrParams.Phis, _phiSuffix);

  // 3.: set magnitudes of all lamlam amplitudes
  setAbsPawianParamsJPCLamLam(upar, theValParams.MagLamLams, theErrParams.MagLamLams, _magSuffix);

  // 4.: set phases of all lam lamlam amplitudes
  setAbsPawianParamsJPCLamLam(upar, theValParams.PhiLamLams, theErrParams.PhiLamLams, _phiSuffix);

  // 5.: set magnitudes of all ls amplitudes
  setAbsPawianParamsLS(upar, theValParams.MagsLS, theErrParams.MagsLS, _magSuffix);

  // 6.: set phases of all ls amplitudes
  setAbsPawianParamsLS(upar, theValParams.PhisLS, theErrParams.PhisLS, _phiSuffix);

   // 7.: set all masses
  setAbsPawianParamsDouble(upar, theValParams.Masses, theErrParams.Masses, _massSuffix);

   // 8.: set all widths
   setAbsPawianParamsDouble(upar, theValParams.Widths, theErrParams.Widths, _widthSuffix);

   // 9.: set all gFactors
  setAbsPawianParamsDouble(upar, theValParams.gFactors, theErrParams.gFactors, _gFactorSuffix);

   // 10.: set all other parameters
  setAbsPawianParamsDouble(upar, theValParams.otherParams, theErrParams.otherParams, _otherParamSuffix);

}



void FitParColBase::printParams(fitParCol& theParams){

  printJPCParams(theParams, theParams.MagsJPC, _magSuffix);
  printJPCParams(theParams, theParams.PhisJPC, _phiSuffix);

  printJPCLSParams(theParams, theParams.Mags, _magSuffix);
  printJPCLSParams(theParams, theParams.Phis, _phiSuffix);
  
  printJPCLamLamParams(theParams, theParams.MagLamLams, _magSuffix);
  printJPCLamLamParams(theParams, theParams.PhiLamLams, _phiSuffix);

  printLSParams(theParams, theParams.MagsLS, _magSuffix);
  printLSParams(theParams, theParams.PhisLS, _phiSuffix);

  printDoubleParams(theParams, theParams.Masses, _massSuffix);
  printDoubleParams(theParams, theParams.Widths, _widthSuffix);
  printDoubleParams(theParams, theParams.gFactors, _gFactorSuffix);
  printDoubleParams(theParams, theParams.otherParams, _otherParamSuffix);
}

void FitParColBase::dumpParams(std::ostream& os, fitParCol& theVals,  fitParCol& theErrs){

  os << std::setprecision(16);
  dumpJPCParams(os, theVals.MagsJPC, theErrs.MagsJPC, _magSuffix);
  dumpJPCParams(os, theVals.PhisJPC, theErrs.PhisJPC, _phiSuffix);
  
  dumpJPCLSParams(os, theVals.Mags, theErrs.Mags,  _magSuffix);
  dumpJPCLSParams(os, theVals.Phis, theErrs.Phis,  _phiSuffix);

  dumpLamLamParams(os, theVals.MagLamLams, theErrs.MagLamLams, _magSuffix);
  dumpLamLamParams(os, theVals.PhiLamLams, theErrs.PhiLamLams, _phiSuffix);

  dumpLSParams(os, theVals.MagsLS, theErrs.MagsLS,  _magSuffix);
  dumpLSParams(os, theVals.PhisLS, theErrs.PhisLS,  _phiSuffix);

  dumpDoubleParams(os, theVals.Masses, theErrs.Masses, _massSuffix);
  dumpDoubleParams(os, theVals.Widths, theErrs.Widths, _widthSuffix);
  dumpDoubleParams(os, theVals.gFactors, theErrs.gFactors, _gFactorSuffix);
  dumpDoubleParams(os, theVals.otherParams, theErrs.otherParams, _otherParamSuffix);

}

void FitParColBase::getFitParamVal(const std::vector<double>& par, fitParCol& theParams){

  unsigned int counter=0;
  getFitParamValJPC(par, theParams.MagsJPC, counter);
  getFitParamValJPC(par, theParams.PhisJPC, counter);
  getFitParamValJPCLS(par, theParams.Mags, counter);
  getFitParamValJPCLS(par, theParams.Phis, counter);
  getFitParamValJPCLamLam(par, theParams.MagLamLams, counter);
  getFitParamValJPCLamLam(par, theParams.PhiLamLams, counter);
  getFitParamValLS(par, theParams.MagsLS, counter);
  getFitParamValLS(par, theParams.PhisLS, counter);
  getFitParamValDouble(par, theParams.Masses, counter);
  getFitParamValDouble(par, theParams.Widths, counter);
  getFitParamValDouble(par, theParams.gFactors, counter);
  getFitParamValDouble(par, theParams.otherParams, counter);

  if(counter!=par.size()){
   Alert << "Parameter vector size does not match the fitParCol parameter count." << endmsg;
  }

}


void FitParColBase::setAbsPawianParamsJPC(std::shared_ptr<AbsPawianParameters> upar, mapStrJPC& startJPCMap, mapStrJPC& errJPCMap, const std::string& suffix){

  mapStrJPC::iterator itJPCMap;
  for (itJPCMap=startJPCMap.begin(); itJPCMap!=startJPCMap.end(); ++itJPCMap){

    std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >& errJPCs= errJPCMap[itJPCMap->first];

    std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::iterator itJPC;
    for (itJPC=itJPCMap->second.begin(); itJPC!=itJPCMap->second.end(); ++itJPC){

      std::shared_ptr<const jpcRes> theJPC=itJPC->first;
      double theStartVal=itJPC->second;
      double theErrVal=errJPCs[theJPC];
      if (theErrVal<=0.) theErrVal=0.1;

      //now fill the fitParameterMap
      std::string jpcStr = theJPC->name() + itJPCMap->first +suffix;

      if (suffix==_phiSuffix){
	 upar->Add(jpcStr, theStartVal, theErrVal);
      }
      else{
	double valMin=0.0;
	double valMax=theStartVal+30.*theErrVal;

	upar->Add(jpcStr, theStartVal, theErrVal, valMin, valMax);
      }
     }

  }

}

void FitParColBase::setAbsPawianParamsLS(std::shared_ptr<AbsPawianParameters> upar, mapStrLS& startMagMap, mapStrLS& errMagMap, const std::string& suffix){

  mapStrLS::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){

    std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& errMags=errMagMap[itMagMap->first];

    std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagMap->second.begin(); itMag!=itMagMap->second.end(); ++itMag){

      std::shared_ptr<const LScomb> theLS=itMag->first;
      double theStartVal=itMag->second;
      double theErrVal=errMags[theLS];
      if (theErrVal<=0.) theErrVal=0.1;

      //now fill the fitParameterMap
      std::string magStr=theLS->name()+itMagMap->first+suffix;

      double valMin=0.;
      double valMax=0.;
      if (suffix==_phiSuffix){
	 //valMin=-4.*M_PI;
	 //valMax=4.*M_PI;
	 upar->Add(magStr, theStartVal, theErrVal);//, valMin, valMax);
      }
      else{
	valMin=theStartVal-6.*theErrVal;
	if (valMin<0.) valMin=0.;
	valMin=0.0;
	valMax=theStartVal+30.*theErrVal;

	upar->Add(magStr, theStartVal, theErrVal, valMin, valMax);
      }
    }

  }

}


void FitParColBase::setAbsPawianParamsJPCLamLam(std::shared_ptr<AbsPawianParameters> upar, mapStrJPCLamLam& startLamLamMagMap, mapStrJPCLamLam& errLamLamMagMap, const std::string& suffix){

  mapStrJPCLamLam::iterator itLamLamMagMap;
  for (itLamLamMagMap=startLamLamMagMap.begin(); itLamLamMagMap!=startLamLamMagMap.end(); ++itLamLamMagMap){

    std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >& errLamLamMags= errLamLamMagMap[itLamLamMagMap->first];

    std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLamMag;
    for (itLamLamMag=itLamLamMagMap->second.begin(); itLamLamMag!=itLamLamMagMap->second.end(); ++itLamLamMag){

      std::shared_ptr<const JPClamlam> theJPCLamLam=itLamLamMag->first;
      double theStartVal=itLamLamMag->second;
      double theErrVal=errLamLamMags[theJPCLamLam];
      if (theErrVal<=0.) theErrVal=0.1;

      //now fill the fitParameterMap
      std::string magStr=theJPCLamLam->name()+itLamLamMagMap->first+suffix;

      double valMin=0.;
      double valMax=0.;
      if (suffix==_phiSuffix){
	 //valMin=-4.*M_PI;
	 //valMax=4.*M_PI;
	 upar->Add(magStr, theStartVal, theErrVal);//, valMin, valMax);
      }
      else{
	valMin=theStartVal-6.*theErrVal;
	if (valMin<0.) valMin=0.;
	valMin = 0.0;
	valMax=theStartVal+30.*theErrVal;

	upar->Add(magStr, theStartVal, theErrVal, valMin, valMax);
      }
    }

  }

}

void FitParColBase::setAbsPawianParamsJPCLS(std::shared_ptr<AbsPawianParameters> upar, mapStrJPCLS& startMagMap, mapStrJPCLS& errMagMap, const std::string& suffix){

  mapStrJPCLS::iterator itMagMap;
  for (itMagMap=startMagMap.begin(); itMagMap!=startMagMap.end(); ++itMagMap){

    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMags=errMagMap[itMagMap->first];

    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=itMagMap->second.begin(); itMag!=itMagMap->second.end(); ++itMag){

      std::shared_ptr<const JPCLS> theJPCLS=itMag->first;
      double theStartVal=itMag->second;
      double theErrVal=errMags[theJPCLS];
      if (theErrVal<=0.) theErrVal=0.1;

      //now fill the fitParameterMap
      std::string magStr=theJPCLS->name()+itMagMap->first+suffix;

      double valMin=0.;
      double valMax=0.;
      if (suffix==_phiSuffix){
	 //valMin=-4.*M_PI;
	 //valMax=4.*M_PI;
	 upar->Add(magStr, theStartVal, theErrVal);//, valMin, valMax);
      }
      else{
	valMin=theStartVal-6.*theErrVal;
	if (valMin<0.) valMin=0.;
	valMin=0.0;
	valMax=theStartVal+30.*theErrVal;

	upar->Add(magStr, theStartVal, theErrVal, valMin, valMax);
      }
    }

  }

}

void FitParColBase::setAbsPawianParamsDouble(std::shared_ptr<AbsPawianParameters> upar, mapStrDouble& startDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix){

  mapStrDouble::iterator it;
  for (it=startDoubleMap.begin(); it!=startDoubleMap.end(); ++it){
    std::string theName =it->first + suffix;
    double theStartVal=it->second;
    double theErrVal=errDoubleMap[it->first];
    if (theErrVal<=0.) theErrVal=0.1;

    double minVal=theStartVal-6.*theErrVal;
    if (minVal<0.) minVal=0.;
    minVal=0;

    double maxVal=theStartVal+30.*theErrVal;

    // for mass and width Parameter
    if (suffix==_massSuffix || suffix==_widthSuffix){
      minVal=theStartVal-5.*theErrVal;
      maxVal=theStartVal+5.*theErrVal;
      if(minVal<0.) minVal=0.;

      upar->Add(theName, theStartVal, theErrVal, minVal, maxVal);
    }
    // for complex fit parameter; phi component; quick workaround
    else if(theName.size()>9 && (theName.compare(theName.size()-8, theName.size(), "PhiOther")==0)){
	minVal=-4.*M_PI;
	maxVal=4.*M_PI;
	upar->Add(theName, theStartVal, theErrVal);//, minVal, maxVal);
    }
    // for parameter where pos and neg values are allowed
    else if(theName.size()>12 && (theName.compare(theName.size()-11, theName.size(), "PosNegOther")==0)){
	minVal = -fabs(theStartVal)-30.*theErrVal;
	maxVal = fabs(theStartVal)+30.*theErrVal;;
	upar->Add(theName, theStartVal, theErrVal);//, minVal, maxVal);
    }
    else if(theName.size()>9 && (theName.compare(theName.size()-8, theName.size(), "PosOther")==0)){
	minVal = 0.;
	maxVal = fabs(theStartVal)+10.*theErrVal;;
	upar->Add(theName, theStartVal, theErrVal);//, minVal, maxVal);
    }
    // for parameter where pos and neg values are allowed
    else if(theName.size()>13 && (theName.compare(theName.size()-12, theName.size(), "Range01Other")==0)){
	minVal = 0.;
	maxVal = 1.;
	upar->Add(theName, theStartVal, theErrVal, minVal, maxVal);
    }
    // channel scaling params
    else if(theName.size()>20 && (theName.compare(theName.size()-19, theName.size(), "channelScalingOther")==0)){
      minVal = 0.;
      //      maxVal = theStartVal+30.*theErrVal;
       upar->Add(theName, theStartVal, theErrVal, minVal, maxVal*10.);
    }
    else{
       upar->Add(theName, theStartVal, theErrVal, minVal, maxVal);
    }
  }

}

void FitParColBase::printJPCParams(fitParCol& theParams, mapStrJPC& jpcMap, const std::string& suffix){

  mapStrJPC::const_iterator it;

  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCMap;

   for (it=jpcMap.begin(); it!=jpcMap.end(); ++it){
     std::string currentJPCName=it->first;

     std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > currentJPCmap=it->second;
     for ( itJPCMap=currentJPCmap.begin(); itJPCMap!=currentJPCmap.end(); ++itJPCMap){
       std::string currentName=itJPCMap->first->name()+currentJPCName+suffix;

       Info << currentName << "\t" << itJPCMap->second << endmsg;
     }
   }
}

void FitParColBase::printLSParams(fitParCol& theParams, mapStrLS& LSMagMap, const std::string& suffix){

  mapStrLS::const_iterator it;

  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::const_iterator itLSMap;

   for (it=LSMagMap.begin(); it!=LSMagMap.end(); ++it){

    std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;
     for ( itLSMap=currentAmp.begin(); itLSMap!=currentAmp.end(); ++itLSMap){
       std::string currentName=itLSMap->first->name()+it->first+suffix;
        Info<< currentName << "\t" << itLSMap->second << endmsg;
     }
   }

}

void FitParColBase::printJPCLamLamParams(fitParCol& theParams, mapStrJPCLamLam& lamLamMap, const std::string& suffix){

  mapStrJPCLamLam::const_iterator it;

  std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLamLamMap;

   for (it=lamLamMap.begin(); it!=lamLamMap.end(); ++it){
     std::string currentAmpName=it->first;

     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;
     for ( itJPCLamLamMap=currentAmp.begin(); itJPCLamLamMap!=currentAmp.end(); ++itJPCLamLamMap){
       std::string currentName=itJPCLamLamMap->first->name()+currentAmpName+suffix;

       Info << currentName << "\t" << itJPCLamLamMap->second << endmsg;
     }
   }
}


void FitParColBase::printJPCLSParams(fitParCol& theParams, mapStrJPCLS& JPCLSMagMap, const std::string& suffix){

  mapStrJPCLS::const_iterator it;

  std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPCLSMap;

   for (it=JPCLSMagMap.begin(); it!=JPCLSMagMap.end(); ++it){

    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;
     for ( itJPCLSMap=currentAmp.begin(); itJPCLSMap!=currentAmp.end(); ++itJPCLSMap){
       std::string currentName=itJPCLSMap->first->name()+it->first+suffix;
        Info<< currentName << "\t" << itJPCLSMap->second << endmsg;
     }
   }

}

void FitParColBase::printDoubleParams(fitParCol& theParams, mapStrDouble& doubleMap, const std::string& suffix){
  mapStrDouble::const_iterator it;

  for (it=doubleMap.begin(); it!=doubleMap.end(); ++it){
    std::string currentName=it->first+suffix;
    Info << currentName << "\t" << it->second;
  }


}

void FitParColBase::getFitParamValJPC(const std::vector<double>& par, mapStrJPC& jpcMap, unsigned int& counter){

  mapStrJPC::iterator it;
  for (it=jpcMap.begin(); it!=jpcMap.end(); ++it){

    std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::iterator itJPC;
    for (itJPC=it->second.begin(); itJPC!=it->second.end(); ++itJPC){
      itJPC->second=par.at(counter);
      counter++;

    }
  }
}

void FitParColBase::getFitParamValLS(const std::vector<double>& par, mapStrLS& lsMap, unsigned int& counter){
  mapStrLS::iterator it;
  for (it=lsMap.begin(); it!=lsMap.end(); ++it){

    std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::iterator itls;
    for (itls=it->second.begin(); itls!=it->second.end(); ++itls){
      itls->second=par.at(counter);
      counter++;

    }
  }
}

void FitParColBase::getFitParamValJPCLamLam(const std::vector<double>& par, mapStrJPCLamLam& lamLamMagMap, unsigned int& counter){

  mapStrJPCLamLam::iterator it;
  for (it=lamLamMagMap.begin(); it!=lamLamMagMap.end(); ++it){

    std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itJPClamlam;
    for (itJPClamlam=it->second.begin(); itJPClamlam!=it->second.end(); ++itJPClamlam){
      itJPClamlam->second=par.at(counter);
      counter++;

    }
  }
}

void FitParColBase::getFitParamValJPCLS(const std::vector<double>& par, mapStrJPCLS& jpclsMap, unsigned int& counter){
  mapStrJPCLS::iterator it;
  for (it=jpclsMap.begin(); it!=jpclsMap.end(); ++it){
    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itjpcls;
    for (itjpcls=it->second.begin(); itjpcls!=it->second.end(); ++itjpcls){
      itjpcls->second=par.at(counter);
      counter++;

    }
  }
}


void FitParColBase::getFitParamValDouble(const std::vector<double>& par, mapStrDouble& doubleMap, unsigned int& counter){
  //  Info << "getFitParamValDouble par[" << counter << "]=\t" << par.at(counter) << endmsg;
  mapStrDouble::iterator it;
  for (it=doubleMap.begin(); it!=doubleMap.end(); ++it){
    it->second=par.at(counter);
    counter++;
  }
}

void  FitParColBase::dumpJPCParams(std::ostream& os, mapStrJPC& valJPCMap, mapStrJPC& errJPCMap, const std::string& suffix){

  mapStrJPC::const_iterator it;

   for (it=valJPCMap.begin(); it!=valJPCMap.end(); ++it){

     std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > currentErrMap=errJPCMap[it->first];

     std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess > currentMap=it->second;

     std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::const_iterator itJPC;
     for ( itJPC=currentMap.begin(); itJPC!=currentMap.end(); ++itJPC){
       std::string currentName=itJPC->first->name()+it->first+suffix;

       os << currentName << "\t" << itJPC->second << "\t" << currentErrMap[itJPC->first] << std::endl;
     }
   }

}

void FitParColBase::dumpLSParams(std::ostream& os, mapStrLS& valLSMap, mapStrLS& errLSMap, const std::string& suffix){

  mapStrLS::const_iterator it;

   for (it=valLSMap.begin(); it!=valLSMap.end(); ++it){

     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentErrMap=errLSMap[it->first];

     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;

     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::const_iterator itLS;
     for ( itLS=currentAmp.begin(); itLS!=currentAmp.end(); ++itLS){
       std::string currentName=itLS->first->name()+it->first+suffix;

       os << currentName << "\t" << itLS->second << "\t" << currentErrMap[itLS->first] << std::endl;
     }
   }
}

void  FitParColBase::dumpLamLamParams(std::ostream& os, mapStrJPCLamLam& valLamLamMap, mapStrJPCLamLam& errLamLamMap, const std::string& suffix){

  mapStrJPCLamLam::const_iterator it;

   for (it=valLamLamMap.begin(); it!=valLamLamMap.end(); ++it){

     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentErrMap=errLamLamMap[it->first];

     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;

     std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::const_iterator itLamLam;
     for ( itLamLam=currentAmp.begin(); itLamLam!=currentAmp.end(); ++itLamLam){
       std::string currentName=itLamLam->first->name()+it->first+suffix;

       os << currentName << "\t" << itLamLam->second << "\t" << currentErrMap[itLamLam->first] << std::endl;
     }
   }

}

void FitParColBase::dumpJPCLSParams(std::ostream& os, mapStrJPCLS& valLSMap, mapStrJPCLS& errLSMap, const std::string& suffix){

  mapStrJPCLS::const_iterator it;

   for (it=valLSMap.begin(); it!=valLSMap.end(); ++it){

     std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentErrMap=errLSMap[it->first];

     std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentAmp=it->second;

     std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::const_iterator itLS;
     for ( itLS=currentAmp.begin(); itLS!=currentAmp.end(); ++itLS){
       std::string currentName=itLS->first->name()+it->first+suffix;

       os << currentName << "\t" << itLS->second << "\t" << currentErrMap[itLS->first] << std::endl;
     }
   }
}

void FitParColBase::dumpDoubleParams(std::ostream& os, mapStrDouble& valDoubleMap, mapStrDouble& errDoubleMap, const std::string& suffix){

  mapStrDouble::const_iterator it;

  for (it=valDoubleMap.begin(); it!=valDoubleMap.end(); ++it){
    std::string currentName=it->first+suffix;
      os << currentName << "\t" << it->second << "\t" << errDoubleMap[it->first] << std::endl;
    }

}
