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

#ifndef _ArgandModuloConverter_H
#define _ArgandModuloConverter_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <memory>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TTree;
class TGraphErrors;
class PhaseSpaceIsobar;


class ArgandModuloConverter {

public:

  // create/copy/destroy:

  ///Constructor 
  ArgandModuloConverter(std::string rootFileInName, double mass1, double mass2, double mass3, double mass4);


  /** Destructor */
  virtual ~ArgandModuloConverter();

  // Getters:

  void convertArgandToModulo();
 
protected:


private:
  TFile* _theTFileOut;
  TFile* _theTFileIn;
  TTree* _dataTree;
  TTree* _fitTree;

  TGraphErrors* _argandUnitsDataGraphErr;
  TGraphErrors* _argandUnitsFitGraphErr;

  TGraphErrors* _moduloDataGraphErr;
  TGraphErrors* _moduloFitGraphErr;

  PhaseSpaceIsobar* _php1;
  PhaseSpaceIsobar* _php2;
  
  float _massVal;
  float _dataVal;
  float _dataErrVal;
  float _fitVal;

  float _moduloDataVal;
  float _moduloDataErrVal;
  float _moduloFitVal;
};

#endif
