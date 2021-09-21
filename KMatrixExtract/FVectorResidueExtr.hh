//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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
//FVectorResidueExtr class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <memory>
#include "math.h" 

#include "KMatrixExtract/TMatrixResidueExtr.hh"
class AbsPawianParameters;
class AbsPhaseSpace;
class PwaCovMatrix;
class pipiScatteringParser;
class FVectorIntensityDynamics;
class FVector;

class FVectorResidueExtr : public TMatrixResidueExtr {

public:

  // create/copy/destroy:

  ///Constructor 
  FVectorResidueExtr(pipiScatteringParser* theParser);

  /** Destructor */
  virtual ~FVectorResidueExtr();

  // Getters:
  virtual void Calculation();
  virtual void CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, std::complex<double>& polePos, std::vector<ResidueProperties>& , std::vector<ResidueProperties>&, std::vector<ResidueProperties>&);
  virtual void fillParams();
  void updateFMatDy(std::shared_ptr<AbsPawianParameters> params);

protected:
  virtual void dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal, std::vector<ResidueProperties> resPropImag, std::vector<ResidueProperties> resPropAv);
  std::shared_ptr<FVectorIntensityDynamics> _fVectorIntensityDynamics;
  std::shared_ptr<FVector> _fVector; 
  std::string _pVecName;
  std::string _extractionMethod;
  std::vector<std::string> _fVecParamNames;

private:
  void init();
};
