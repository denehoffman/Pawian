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

#include <fstream>

#include "FitParams/StreamFitParColBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"

StreamFitParmsBase::StreamFitParmsBase(std::string& filePath) :
  AbsFitParColStreamer(filePath)
{
  _paramVal = GlobalEnv::instance()->DefaultParamVal();
  _paramErr = GlobalEnv::instance()->DefaultParamErr();
  fillParams();
}

StreamFitParmsBase::~StreamFitParmsBase(){;}

void StreamFitParmsBase::fillParams(){
  const std::string magSuffix="Mag";
  const std::string phiSuffix="Phi";
  const std::string massSuffix="Mass";
  const std::string widthSuffix="Width";
  const std::string gFactorSuffix="gFactor";
  const std::string otherSuffix="Other";

  fillJPCAmps(_paramVal.MagsJPC, _paramErr.MagsJPC, magSuffix);
  fillJPCAmps(_paramVal.PhisJPC, _paramErr.PhisJPC, phiSuffix);
  fillLSAmps(_paramVal.MagsLS, _paramErr.MagsLS, magSuffix);
  fillLSAmps(_paramVal.PhisLS, _paramErr.PhisLS, phiSuffix);
  fillJPCLamLamAmps(_paramVal.MagLamLams, _paramErr.MagLamLams, magSuffix);
  fillJPCLamLamAmps(_paramVal.PhiLamLams, _paramErr.PhiLamLams, phiSuffix);
  fillJPCLSAmps(_paramVal.Mags, _paramErr.Mags, magSuffix);
  fillJPCLSAmps(_paramVal.Phis, _paramErr.Phis, phiSuffix);
  fillDoubles(_paramVal.Masses, _paramErr.Masses, massSuffix);
  fillDoubles(_paramVal.Widths, _paramErr.Widths, widthSuffix);
  fillDoubles(_paramVal.gFactors, _paramErr.gFactors, gFactorSuffix);
  fillDoubles(_paramVal.otherParams, _paramErr.otherParams, otherSuffix);
}

void StreamFitParmsBase::fillJPCAmps(mapStrJPC& valMap, mapStrJPC& errMap,
				     const std::string& suffix){
  mapStrJPC::iterator itJPCMap;
  for( itJPCMap=valMap.begin(); itJPCMap!=valMap.end(); ++itJPCMap){
    std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::iterator itJPC;
    for ( itJPC=itJPCMap->second.begin(); itJPC!=itJPCMap->second.end();  ++itJPC){
      std::string theKey=itJPC->first->name()+itJPCMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;

      stringPairIter=_stringPairMap.find(theKey);

      if ( stringPairIter != _stringPairMap.end() ){
	Info << "fillJPCs: key\t" << theKey << "\tfound" << endmsg;
	double val=stringPairIter->second.first;
	Info << "replace val by " << val << endmsg;
	double err=stringPairIter->second.second;
	Info << "replace err by " << err << endmsg;

	valMap[itJPCMap->first][itJPC->first] = val;
	errMap[itJPCMap->first][itJPC->first] = err;

      }
    }
  }
}

void StreamFitParmsBase::fillLSAmps(mapStrLS& valMap, mapStrLS& errMap, const std::string& suffix){

  mapStrLS::iterator itLSMap;
  for( itLSMap=valMap.begin(); itLSMap!=valMap.end(); ++itLSMap){

    std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::iterator itLS;

    for ( itLS=itLSMap->second.begin(); itLS!=itLSMap->second.end();  ++itLS){


      std::string theKey=itLS->first->name()+itLSMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;

      stringPairIter=_stringPairMap.find(theKey);

      if ( stringPairIter != _stringPairMap.end() ){
	Info << "fillLSAmps: key\t" << theKey << "\tfound" << endmsg;
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

void StreamFitParmsBase::fillJPCLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap, const std::string& suffix){

  mapStrJPCLamLam::iterator itLamLamMap;
  for( itLamLamMap=valMap.begin(); itLamLamMap!=valMap.end(); ++itLamLamMap){

    std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLam;

    for ( itLamLam=itLamLamMap->second.begin(); itLamLam!=itLamLamMap->second.end();  ++itLamLam){


      std::string theKey=itLamLam->first->name()+itLamLamMap->first+suffix;
      Info << "fillJPCLamLamAmps: theKey=\t" << theKey << endmsg;
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

void StreamFitParmsBase::fillJPCLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap, const std::string& suffix){

  mapStrJPCLS::iterator itLSMap;
  for( itLSMap=valMap.begin(); itLSMap!=valMap.end(); ++itLSMap){

    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itLS;

    for ( itLS=itLSMap->second.begin(); itLS!=itLSMap->second.end();  ++itLS){


      std::string theKey=itLS->first->name()+itLSMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;

      stringPairIter=_stringPairMap.find(theKey);

      if ( stringPairIter != _stringPairMap.end() ){
	Info << "fillJPCLSAmps: key\t" << theKey << "\tfound" << endmsg;
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

void StreamFitParmsBase::fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap, const std::string& suffix){

  mapStrDouble::iterator itMap;
  for( itMap=valMap.begin(); itMap!=valMap.end(); ++itMap){

    std::string theKey=itMap->first+suffix;
    Info << "fillDoubles: theKey=\t" << theKey << endmsg;
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

void StreamFitParmsBase::fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index){

  StringPairMap::const_iterator stringPairIter=_stringPairMap.find(suffix);
  if ( stringPairIter != _stringPairMap.end() ){
    theValMap[index]=stringPairIter->second.first;
    theErrMap[index]=stringPairIter->second.second;
  }

}
