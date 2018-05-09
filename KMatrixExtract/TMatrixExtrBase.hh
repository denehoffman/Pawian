//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//TMatrixExtrBase class definition file. -*- C++ -*-
// Copyright 2018 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <memory>

class KMatrixParser;
class AbsPhaseSpace;
class TMatrixRel;
class KMatrixRel;
class KPole;
class ParticleTable;
class KMatrixParser;
class AbsPawianParameters;

class TMatrixExtrBase {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixExtrBase(std::string pathToConfigParser, std::string pathToFitParams);


  /** Destructor */
  virtual ~TMatrixExtrBase();

  // Getters:

protected:
  std::shared_ptr<KMatrixParser> _kMatrixParser;
  std::vector< std::string> _gFactorNames;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  std::shared_ptr<TMatrixRel> _tMatr;
  std::shared_ptr<KMatrixRel> _kMatr;
  ParticleTable* _particleTable;
  
  std::shared_ptr<AbsPawianParameters> _params;
  std::string _pathToFitParams;
  int _orbitalL;
  std::vector<double> _signs;

private:
  void init();
};


