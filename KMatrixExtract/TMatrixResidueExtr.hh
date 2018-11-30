//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//                 Xiaoshuai Qin (xqin@ep1.rub.de)  
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

#include "KMatrixExtract/TMatrixExtrBase.hh"
class AbsPawianParameters;
class AbsPhaseSpace;
class PwaCovMatrix;
class TMatrixExtrFit;

struct ResidueProperties {
  double absR;
  double theta;
  double gammai;
  double errAbsR;
  double errTheta;
  double errGammai;
};

class TMatrixResidueExtr : public TMatrixExtrBase {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixResidueExtr(std::string pathToConfigParser, std::string pathToFitParams, std::string sheet, std::string pathToSerialzationFile, std::complex<double> energyBorderMin, std::complex<double> energyBorderMax, std::complex<double> energyStartParams);


  /** Destructor */
  virtual ~TMatrixResidueExtr();

  // Getters:
  void Calculation();
  void CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, std::complex<double>& polePos, std::vector<ResidueProperties>& , std::vector<ResidueProperties>&, std::vector<ResidueProperties>&);
  //  std::complex<double>  CalcMassWidth();
  std::complex<double> CalcMassWidth(std::shared_ptr<AbsPawianParameters> currentParameters);
  bool GetCovMatrix();

protected:


private:
  void dumpResult(std::complex<double> polePos, std::vector<ResidueProperties> resPropReal, std::vector<ResidueProperties> resPropImag, std::vector<ResidueProperties> resPropAv);

  std::string _pathToSerialzationFile;
  std::complex<double> _energyMin;
  std::complex<double> _energyMax;
  std::complex<double> _energyStart;
  std::shared_ptr<TMatrixExtrFit> _tMatFit;
  std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
};
