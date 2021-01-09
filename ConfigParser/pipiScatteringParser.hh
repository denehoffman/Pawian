//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pipiScatteringParser class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf, 2020 Malte Albrecht

#pragma once

#include "ConfigParser/ParserBase.hh"
// Boost headers go here



class pipiScatteringParser : public ParserBase 
{

  public:

  pipiScatteringParser(int argc,char **argv);
  virtual ~pipiScatteringParser(){;}

  const std::string pathToKMatrixCompareFile() const {return _pathToKMatrixCompareFile;}
  const std::string baseNameFVector() const {return _baseNameFVector;}
  const std::string baseNameFVectorCompare() const {return _baseNameFVectorCompare;}
  const double minImagMass() const {return _minImagMass;}
  const double maxImagMass() const {return _maxImagMass;}
  const double minRealMass() const {return _minRealMass;}
  const double maxRealMass() const {return _maxRealMass;}
  const double startRealMass() const {return _startRealMass;}
  const double startImagMass() const {return _startImagMass;}
  const int noOfStepsForScan() {return _noOfStepsForScan;}
  const std::string sheet() const {return _sheet;}
  const std::string prodMomParamsStr() const {return _prodMomParamsStr;}

protected:
  virtual bool parseCommandLine(int argc,char **argv); 
  std::string _pathToKMatrixCompareFile;
  std::string _baseNameFVector;
  std::string _baseNameFVectorCompare;
  double _minImagMass;
  double _maxImagMass;
  double _minRealMass;
  double _maxRealMass;
  double _startRealMass;
  double _startImagMass;
  int _noOfStepsForScan;
  std::string _sheet;
  std::string _prodMomParamsStr;
};


