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

#include "PwaUtils/AbsFitParamStreamer.hh"
#include "ErrLogger/ErrLogger.hh"

AbsFitParamStreamer::AbsFitParamStreamer(std::string& filePath)
{

  std::ifstream ifs(filePath.c_str());
  if(!ifs.good()) 
    {
      Alert << "Error: could not parse " << filePath << endmsg;
      exit(1);
    }


  std::string strTmp, strState;
  bool bError=false;
  double dMag=0.0, dPhi=0.0;
  
  while(!ifs.eof())
    {
      ifs >> strState;
      if (ifs.eof()) break;
      ifs >> strTmp;
      dMag = atof(strTmp.c_str());
      if (ifs.eof()) break;
      ifs >> strTmp;
      dPhi = atof(strTmp.c_str());
      if (ifs.eof()) break;
      if (!bError) 
	{
	  DebugMsg << "Start Parameter: " << strState << " Mag: " << dMag << " Phi: " << dPhi << " set.\n" << endmsg;
	  StringParameterMap(strState,dMag,dPhi);
	}
    }


}

AbsFitParamStreamer::~AbsFitParamStreamer(){;}


void AbsFitParamStreamer::StringParameterMap(const std::string &theName, const double firstVal, const double secondVal)
{
  std::pair <double,double> thePair=make_pair(firstVal, secondVal);
  _stringPairMap[theName]=thePair;
}

void AbsFitParamStreamer::fillParamMapAmps(std::vector< std::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< std::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& toFill){

  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;

  for ( itJPCLS=theJPCLSs.begin(); itJPCLS!=theJPCLSs.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string theKey=(*itJPCLS)->name()+suffix;

    stringPairIter=_stringPairMap.find(theKey);

    if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[theKey];  
    toFill[*itJPCLS]=theValPair;
    DebugMsg << theKey << " filled with values: "  << theValPair.first <<"\t" << theValPair.second  << endmsg;
    }
    else{
      Warning << theKey <<" not available in fitVal input file!!!!"  << endmsg;
      std::pair <double,double> theValPair=make_pair(-10000., -10000.);
      toFill[*itJPCLS]=theValPair;      
    }

  }
}

void AbsFitParamStreamer::fillParamMapMass(std::string& name, pair<double, double>& toFill){

  StringPairMap::const_iterator stringPairIter;
  stringPairIter=_stringPairMap.find(name);

    if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[name];  
    toFill=theValPair;
    DebugMsg << name << " filled with values: "  << theValPair.first <<"\t" << theValPair.second  << endmsg;
    }
    else{
      Warning << name <<" not available in fitVal input file!!!!"  << endmsg;
      std::pair <double,double> theValPair=make_pair(-10000., -10000.);
      toFill=theValPair;
    }

}

void AbsFitParamStreamer::fillParamFlatte(std::string& name, double& toFill){
  StringPairMap::const_iterator stringPairIter;
  stringPairIter=_stringPairMap.find(name);
  if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[name];
    toFill=theValPair.first;
    DebugMsg << name << " filled with value: "  << theValPair.first <<"\t" << endmsg;
  }
  else{
    Warning << name <<" not available in fitVal input file!!!!"  << endmsg;
    toFill=-10000.;
  }  
}
