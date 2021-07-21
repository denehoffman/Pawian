//************************************************************************//
//									  //
//  Copyright 2021 Meike Kuessner (mkussner@ep1.rub.de)                   //
//                 Bertram Kopf (bertram@ep1.rub.de)			  //
//                 Xiaoshuai Qin (xqin@ep1.rub.de)                        //
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
//TMatrixResiduePathExtr class definition file. -*- C++ -*-
// Copyright 2021 Meike Kuessner, Bertram Kopf

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



class TMatrixResiduePathExtr : public TMatrixResidueExtr {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixResiduePathExtr(pipiScatteringParser* theParser);

  /** Destructor */
  virtual ~TMatrixResiduePathExtr();

  // Getters:
  virtual void CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, std::complex<double>& polePos, std::vector<ResidueProperties>& , std::vector<ResidueProperties>&, std::vector<ResidueProperties>&);

protected:

  
private:

};
