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

// ParamDepGFacToFixGFacsWidthMass class definition file. -*- C++ -*-
// Copyright 2014 Julian Pychy

#include "FitParams/ParamDep.hh"

#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <map>

class AbsPawianParameters; 
class KMatrixParser;

struct FixGFacsWidthMassData
{
  std::string poleMass;
  std::string gFacReleaseName;
  std::vector<std::string> gFacScaleNames;
  std::vector<double> massCorPol;
  std::vector<double> widthCorPol;
  std::map<std::string, double> gFacSeparateScales;  

  void print(std::ostream& os) const {
    os <<"pole: " << poleMass
       <<"\t released g-factor: " << gFacReleaseName << std::endl;

    os <<"g factors to be scaled:";
    std::vector<std::string>::const_iterator itStr;

    for(itStr=gFacScaleNames.begin(); itStr != gFacScaleNames.end(); ++itStr){
      os << "\t" << (*itStr);
    } 

    os << "\nparameter for pole mass correction: ";
    std::vector<double>::const_iterator itDouble;
    for(itDouble=massCorPol.begin(); itDouble!=massCorPol.end(); ++itDouble){    
      os<<"\t" << (*itDouble);
    }

    os << "\nparameter for pole width correction: ";
    for(itDouble=widthCorPol.begin(); itDouble!=widthCorPol.end(); ++itDouble){    
      os<<"\t" << (*itDouble);
    }

    os << "\ng factor with separate scaling to released g-factor: ";
    std::map<std::string, double>::const_iterator itMap;
    for(itMap=gFacSeparateScales.begin(); itMap!=gFacSeparateScales.end(); ++itMap){    
      os<<"\t" << itMap->first << "\t factor to " << gFacReleaseName << "\t" << itMap->second << std::endl;
    }

    os<<"\n" << std::endl;
  }
};


class ParamDepGFacToFixGFacsWidthMass : public ParamDep
{
public:
  ParamDepGFacToFixGFacsWidthMass(std::istringstream& configLine, std::shared_ptr<AbsPawianParameters> params);
  ~ParamDepGFacToFixGFacsWidthMass();

  virtual void FillDerived(std::istringstream& configLine);
  virtual void Apply(std::shared_ptr<AbsPawianParameters> params);

private:
  double _targetFullWidth;
  std::shared_ptr<KMatrixParser> _kMatrixParser;
  std::vector<std::string> _poleNames;
  std::map<std::string, FixGFacsWidthMassData> _dataMap;
  std::map<std::string, double> _gFactorStartMap;
};



