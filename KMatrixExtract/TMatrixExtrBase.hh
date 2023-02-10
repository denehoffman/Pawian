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

#include "KMatrixExtract/TMatrixGeneralBase.hh"

class KMatrixParser;
class AbsPhaseSpace;
class TMatrixRel;
class KMatrixRel;
class KPole;
class ParticleTable;
class KMatrixParser;
class AbsPawianParameters;
class TMatrixDynamics;

class TMatrixExtrBase : public TMatrixGeneralBase{

public:

  // create/copy/destroy:

  ///Constructor
  TMatrixExtrBase(pipiScatteringParser* theParser);

  /** Destructor */
  virtual ~TMatrixExtrBase();
  virtual double calcTMatrix(double eReal, double eImag);
  void updateTMatDy(std::shared_ptr<AbsPawianParameters> params);
  std::shared_ptr<TMatrixRel> getNewTMat();
  virtual void SetParamValue(const std::string & paramName, double paramVal);
  void setSheet(std::string newSheet);
  // Getters:

protected:
  std::string _sheet;
  std::vector<double> _signs;

private:
// void init();
};


