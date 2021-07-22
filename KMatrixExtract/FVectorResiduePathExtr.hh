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
//FVectorResiduePathExtr class definition file. -*- C++ -*-
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

#include "KMatrixExtract/FVectorResidueExtr.hh"
class AbsPawianParameters;
class AbsPhaseSpace;
class PwaCovMatrix;
class pipiScatteringParser;
class FVectorIntensityDynamics;
class FVector;

class FVectorResiduePathExtr : public FVectorResidueExtr {

public:

  // create/copy/destroy:

  ///Constructor 
  FVectorResiduePathExtr(pipiScatteringParser* theParser);

  /** Destructor */
  virtual ~FVectorResiduePathExtr();

  // Getters:
  virtual void CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, std::complex<double>& polePos, std::vector<ResidueProperties>& , std::vector<ResidueProperties>&, std::vector<ResidueProperties>&);

protected:

private:
};
