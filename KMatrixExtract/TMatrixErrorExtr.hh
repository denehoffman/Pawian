//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
//                 Xiaoshuai Qin (xqin@ep1.rub.de)             //
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
//TMatrixErrorExtr class definition file. -*- C++ -*-
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
class PwaCovMatrix;
class TMatrixExtrFit;

class TMatrixErrorExtr : public TMatrixExtrBase , public std::enable_shared_from_this<TMatrixErrorExtr> {

public:

  // create/copy/destroy:

  ///Constructor 
  TMatrixErrorExtr(pipiScatteringParser* theParser);

  /** Destructor */
  virtual ~TMatrixErrorExtr();

  // Getters:
  virtual bool GetCovMatrix();
  virtual void CalcOriginal();
  virtual void CalcWithErrrors();
  virtual void Calculation();
  virtual void printErrors();
  virtual std::complex<double> CalcMassWidth(std::shared_ptr<AbsPawianParameters> currentParameters);
  std::complex<double> GetResult() {return _result;}
  std::complex<double> GetError() {return _error;}
  std::complex<double> energyMin() {return _energyMin;}
  std::complex<double> energyMax() {return _energyMax;}
  std::complex<double> energyStart() {return _energyStart;}

protected:
  static bool cmp(std::pair<std::string,double> const & a, std::pair<std::string,double> const & b) 
  { 
	return a.second != b.second?  abs(a.second) > abs(b.second) : a.first < b.first;
  }
  std::string _pathToSerialzationFile;
  std::complex<double> _energyMin;
  std::complex<double> _energyMax;
  std::complex<double> _energyStart;
  std::complex<double> _result;
  std::complex<double> _error;
  std::shared_ptr<PwaCovMatrix> _thePwaCovMatrix;
  std::map< std::string, std::complex<double> > _derivatives;
  std::vector< std::pair<std::string, double> > _realDerivatives;
  std::vector< std::pair<std::string, double> > _imagDerivatives;
  std::vector< std::pair<std::string, double> > _realError;
  std::vector< std::pair<std::string, double> > _imagError;
  bool _calcWithErrors;

private:

};
