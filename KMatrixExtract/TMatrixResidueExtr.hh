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
//TMatrixResidueExtr class definition file. -*- C++ -*-
// Copyright 2018 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <memory>
#include "math.h" 

class AbsPawianParameters;
class AbsPhaseSpace;
class PwaCovMatrix;
class TMatrixExtrFit;

class TMatrixResidueExtr {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixResidueExtr(std::string pathToConfigParser, std::string pathToFitParams, std::string sheet, std::string pathToSerialzationFile, std::complex<double> energyBorderMin, std::complex<double> energyBorderMax, std::complex<double> energyStartParams);


  /** Destructor */
  virtual ~TMatrixResidueExtr();

  // Getters:
  std::complex<double>  CalcResidue();
  void CalcResidueAll();
  std::complex<double>  CalcMassWidth();
  double  calcPartialWidth(double gFac, std::complex<double> poleMass, std::shared_ptr<AbsPhaseSpace> php);
  void GetCovMatrix();

protected:


private:
  std::string _pathToSerialzationFile;
  std::complex<double> _energyMin;
  std::complex<double> _energyMax;
  std::complex<double> _energyStart;
  std::shared_ptr<TMatrixExtrFit> _tMatFit;
  std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
};
