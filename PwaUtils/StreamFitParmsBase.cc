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

#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"

StreamFitParmsBase::StreamFitParmsBase(std::string& filePath) :
  AbsFitParamStreamer(filePath)
{
  _paramVal = GlobalEnv::instance()->DefaultParamVal();
  _paramErr = GlobalEnv::instance()->DefaultParamErr();
  fillParams();
}

StreamFitParmsBase::~StreamFitParmsBase(){;}

void StreamFitParmsBase::fillParams(){
  const std::string isoSuffix="Iso";
  const std::string magSuffix="Mag";
  const std::string phiSuffix="Phi";
  const std::string massSuffix="Mass";
  const std::string widthSuffix="Width";
  const std::string gFactorSuffix="gFactor";
  const std::string otherSuffix="Other";

  fillJPCIsos(_paramVal.Isos, _paramErr.Isos, isoSuffix);
  fillLamLamAmps(_paramVal.MagLamLams, _paramErr.MagLamLams, magSuffix);
  fillLamLamAmps(_paramVal.PhiLamLams, _paramErr.PhiLamLams, phiSuffix);
  fillLSAmps(_paramVal.Mags, _paramErr.Mags, magSuffix);
  fillLSAmps(_paramVal.Phis, _paramErr.Phis, phiSuffix);
  fillDoubles(_paramVal.Masses, _paramErr.Masses, massSuffix);
  fillDoubles(_paramVal.Widths, _paramErr.Widths, widthSuffix);
  fillDoubles(_paramVal.gFactors, _paramErr.gFactors, gFactorSuffix);
  fillDoubles(_paramVal.otherParams, _paramErr.otherParams, otherSuffix);
}

void StreamFitParmsBase::fillJPCIsos(mapStrJPC& valMap, mapStrJPC& errMap,
				     const std::string& suffix){
  mapStrJPC::iterator itIsoMap;
  for( itIsoMap=valMap.begin(); itIsoMap!=valMap.end(); ++itIsoMap){
    std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >::iterator itIso;
    for ( itIso=itIsoMap->second.begin(); itIso!=itIsoMap->second.end();  ++itIso){
      std::string theKey=itIso->first->name()+itIsoMap->first+suffix;
      Info << "theKey=\t" << theKey << endmsg;
      StringPairMap::const_iterator stringPairIter;

      stringPairIter=_stringPairMap.find(theKey);

      if ( stringPairIter != _stringPairMap.end() ){
	Info << "key\t" << theKey << "\tfound" << endmsg;
	double val=stringPairIter->second.first;
	Info << "replace val by " << val << endmsg;
	double err=stringPairIter->second.second;
	Info << "replace err by " << err << endmsg;

	valMap[itIsoMap->first][itIso->first] = val;
	errMap[itIsoMap->first][itIso->first] = err;

      }
    }
  }
}

void StreamFitParmsBase::fillLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap, const std::string& suffix){

  mapStrJPCLamLam::iterator itLamLamMap;
  for( itLamLamMap=valMap.begin(); itLamLamMap!=valMap.end(); ++itLamLamMap){

    std::map< std::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itLamLam;

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

void StreamFitParmsBase::fillLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap, const std::string& suffix){

  mapStrJPCLS::iterator itLSMap;
  for( itLSMap=valMap.begin(); itLSMap!=valMap.end(); ++itLSMap){

    std::map< std::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itLS;

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

void StreamFitParmsBase::fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap, const std::string& suffix){

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

void StreamFitParmsBase::fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index){

  StringPairMap::const_iterator stringPairIter=_stringPairMap.find(suffix);
  if ( stringPairIter != _stringPairMap.end() ){
    theValMap[index]=stringPairIter->second.first;
    theErrMap[index]=stringPairIter->second.second;
  }

}
