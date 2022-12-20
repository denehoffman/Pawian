//************************************************************************//
//                                                                        //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 Julian Pychy (julian@ep1.rub.de)                       //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// Copyright 2014 Julian Pychy

#include <iostream>
#include <sstream>

#include "FitParams/ParamDepGFacToFixGFacsWidthMass.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"

ParamDepGFacToFixGFacsWidthMass::ParamDepGFacToFixGFacsWidthMass(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params){
  std::vector<std::string> targetParameterVec;
  std::string pathToKMatrixConfig;
  configLine >> pathToKMatrixConfig;
  std::string kMatStorePath=getenv("KMAT_DIR");
  std::string completePathToKMatrixConfig=kMatStorePath+"/"+pathToKMatrixConfig;
  std::ifstream ifs(completePathToKMatrixConfig);
  if(!ifs.good()){
    WarningMsg << "file " << completePathToKMatrixConfig << " does not exist. Try now to search in absolute path" << endmsg;
    completePathToKMatrixConfig=pathToKMatrixConfig;
    ifs=std::ifstream(completePathToKMatrixConfig);
    if(!ifs.good()){
      Alert << "file: " << completePathToKMatrixConfig << " does not exist!!!!" << endmsg;
      exit(1);
    }
  }
  

  _kMatrixParser=std::shared_ptr<KMatrixParser>(new KMatrixParser(completePathToKMatrixConfig));
  _gFactorStartMap=_kMatrixParser->gFactorMap();

  const std::vector<std::string> gFactorFixPoles= _kMatrixParser->gFactorFixPoles();
  std::vector<std::string>::const_iterator it;

  for(it=gFactorFixPoles.begin(); it!=gFactorFixPoles.end(); ++it){
    InfoMsg << (*it) <<endmsg;
    FixGFacsWidthMassData newData;
    std::string poleMassName=(*it)+"Mass";
    newData.poleMass=poleMassName;
    targetParameterVec.push_back(*it);
    _dataMap[*it]=newData;
  }

  const std::vector<std::string> gFactorFixReleaseGFactor = _kMatrixParser->gFactorFixReleaseGFactor();
  for(it=gFactorFixReleaseGFactor.begin(); it!=gFactorFixReleaseGFactor.end(); ++it){
    std::istringstream currentStream(*it);
    std::string currentPoleName;
    currentStream >> currentPoleName;
    std::string currentReleasgFactorName;
    currentStream >> currentReleasgFactorName;
    _dataMap.at(currentPoleName).gFacReleaseName=currentReleasgFactorName;       
  }

  const std::vector<std::string> gFactorFixKeepGRatio = _kMatrixParser->gFactorFixKeepGRatio();
  for(it=gFactorFixKeepGRatio.begin(); it!=gFactorFixKeepGRatio.end(); ++it){
     std::istringstream currentStream(*it);
    std::string currentPoleName;
    currentStream >> currentPoleName;
    while (!currentStream.eof()) {
      std::string currentToBeScaledgFactorName;
      currentStream >> currentToBeScaledgFactorName;
      _dataMap.at(currentPoleName).gFacScaleNames.push_back(currentToBeScaledgFactorName);
      targetParameterVec.push_back(currentToBeScaledgFactorName);
    }       
  }

  const std::vector<std::string> gFactorFixMassPol = _kMatrixParser->gFactorFixMassPol();
  for(it=gFactorFixMassPol.begin(); it!=gFactorFixMassPol.end(); ++it){
    std::istringstream currentStream(*it);
    std::string currentPoleName;
    currentStream >> currentPoleName;

    _dataMap.at(currentPoleName).massCorPol.resize(4);
    std::string p0Str;
    currentStream >> p0Str;
    _dataMap.at(currentPoleName).massCorPol[0] = atof(p0Str.c_str());

    std::string p1Str;
    currentStream >> p1Str;
    _dataMap.at(currentPoleName).massCorPol[1] = atof(p1Str.c_str());

    std::string p2Str;
    currentStream >> p2Str;
    _dataMap.at(currentPoleName).massCorPol[2] = atof(p2Str.c_str());

    std::string p3Str;
    currentStream >> p3Str;
    _dataMap.at(currentPoleName).massCorPol[3] = atof(p3Str.c_str());
  }


  const std::vector<std::string> gFactorFixWidthPol = _kMatrixParser->gFactorFixWidthPol();
  for(it=gFactorFixWidthPol.begin(); it!=gFactorFixWidthPol.end(); ++it){
    std::istringstream currentStream(*it);
    std::string currentPoleName;
    currentStream >> currentPoleName;

    _dataMap.at(currentPoleName).widthCorPol.resize(4);
    std::string p0Str;
    currentStream >> p0Str;
    _dataMap.at(currentPoleName).widthCorPol[0] = atof(p0Str.c_str());

    std::string p1Str;
    currentStream >> p1Str;
    _dataMap.at(currentPoleName).widthCorPol[1] = atof(p1Str.c_str());

    std::string p2Str;
    currentStream >> p2Str;
    _dataMap.at(currentPoleName).widthCorPol[2] = atof(p2Str.c_str());

    std::string p3Str;
    currentStream >> p3Str;
    _dataMap.at(currentPoleName).widthCorPol[3] = atof(p3Str.c_str());
  }

  const std::vector<std::string> gFactorFixSeparateScale = _kMatrixParser->gFactorFixSeparateScale();
  for(it=gFactorFixSeparateScale.begin(); it!=gFactorFixSeparateScale.end(); ++it){
    std::istringstream currentStream(*it);
    std::string currentPoleName;
    currentStream >> currentPoleName;

    std::string currentgFacNameToBeScaled;
    currentStream >> currentgFacNameToBeScaled;

    std::string scaleStr;
    currentStream >> scaleStr;

    _dataMap.at(currentPoleName).gFacSeparateScales[currentgFacNameToBeScaled] = atof(scaleStr.c_str());
    targetParameterVec.push_back(currentgFacNameToBeScaled);

  }

  Fill(targetParameterVec, params);

  InfoMsg << "*****ParamDepGFacToFixGFacsWidthMass*****" << endmsg; 
  std::map<std::string, FixGFacsWidthMassData>::iterator itMap;
  for(itMap=_dataMap.begin(); itMap!=_dataMap.end(); ++itMap) itMap->second.print(std::cout);  
}

ParamDepGFacToFixGFacsWidthMass::~ParamDepGFacToFixGFacsWidthMass(){
}

void ParamDepGFacToFixGFacsWidthMass::FillDerived(std::istringstream& configLine){
}

void ParamDepGFacToFixGFacsWidthMass::Apply(std::shared_ptr<AbsPawianParameters> params){

  std::map<std::string, FixGFacsWidthMassData>::iterator itMap;
  for(itMap=_dataMap.begin(); itMap!=_dataMap.end(); ++itMap){
    std::string currentPoleMassName= itMap->second.poleMass;

    double currengFacVal=params->Value(itMap->second.gFacReleaseName);

    double newKpoleMass=itMap->second.massCorPol.at(0)
      +itMap->second.massCorPol.at(1)*currengFacVal
      +itMap->second.massCorPol.at(2)*currengFacVal*currengFacVal
      +itMap->second.massCorPol.at(3)*currengFacVal*currengFacVal*currengFacVal;

    params->SetValue(currentPoleMassName, newKpoleMass);

    double newgFacorScaler=itMap->second.widthCorPol.at(0)
      +itMap->second.widthCorPol.at(1)*currengFacVal
      +itMap->second.widthCorPol.at(2)*currengFacVal*currengFacVal
      +itMap->second.widthCorPol.at(3)*currengFacVal*currengFacVal*currengFacVal;

    std::vector<std::string>::iterator itStr;
    for(itStr=itMap->second.gFacScaleNames.begin(); itStr!=itMap->second.gFacScaleNames.end(); ++itStr){
      //      double newVal= params->Value(*itStr)*newgFacorScaler;
      double newVal= _gFactorStartMap.at(*itStr)*newgFacorScaler;
      params->SetValue((*itStr), newVal);
    } 
   
    std::map<std::string, double>::iterator itStrD;
    for (itStrD=itMap->second.gFacSeparateScales.begin(); itStrD!=itMap->second.gFacSeparateScales.end(); ++itStrD){
      std::string currentgFacName=itStrD->first;
      double newVal=currengFacVal*itStrD->second;
      params->SetValue(currentgFacName, newVal);
    } 
 
  }
}




